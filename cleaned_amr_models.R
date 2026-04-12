
# AMR resistance analysis 
# Refactored from the original analysis script into reusable, GitHub-friendly functions.
# Author: Dana Itani
# Notes:
# - Keeps core modelling logic used in the original script.
# - Removes repeated helper functions and species-specific duplication.
# - Keeps KPN an example same codes has been used for all pathogens ONLY changing abx 

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(readxl)
  library(writexl)
  library(car)
})

# ============================================================
# 1) Label helpers
# ============================================================

drug_label <- function(code) {
  labmap <- c(
    amk = "Amikacin",
    amp = "Ampicillin",
    atm = "Aztreonam",
    caz = "Ceftazidime",
    cip = "Ciprofloxacin",
    ctx = "Cefotaxime",
    etp = "Ertapenem",
    fep = "Cefepime",
    fof = "Fosfomycin",
    gen = "Gentamicin",
    ipm = "Imipenem",
    mem = "Meropenem",
    mero = "Meropenem",
    sxt = "Trimethoprim–Sulfamethoxazole",
    tgc = "Tigecycline",
    tzp = "Piperacillin–Tazobactam"
  )

  code <- tolower(as.character(code))
  out <- unname(labmap[code])
  ifelse(is.na(out), toupper(code), out)
}

bug_label <- function(bug) {
  bug <- toupper(as.character(bug))
  recode(
    bug,
    ECOL = "E. coli",
    KPN  = "K. pneumoniae",
    PSA  = "P. aeruginosa",
    ACN  = "Acinetobacter spp.",
    .default = bug
  )
}

# ============================================================
# 2) AST parsing helpers
# ============================================================

clean_ast <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == ""] <- NA_character_
  x
}

ast_bucket <- function(x) {
  x <- clean_ast(x)

  # Drop numeric-only values unless explicitly mapped elsewhere.
  has_letter <- grepl("[A-Z]", x)
  x[!has_letter] <- NA_character_
  x <- gsub("\\s+", " ", x)

  is_r <- grepl("^R", x) | grepl("RESIST", x) | grepl("NON-?SUSC", x) | grepl("NONSUSC", x)
  is_i <- grepl("^I", x) | grepl("INTER", x) | grepl("^SDD", x) | grepl("DOSE.?DEPENDENT", x)
  is_s <- grepl("^S", x) | grepl("SUSC", x) | grepl("SENSIT", x)

  out <- rep(NA_character_, length(x))
  out[is_s] <- "S"
  out[is_i] <- "I"
  out[is_r] <- "R"
  out
}

safe_relevel <- function(x, ref) {
  x <- factor(x)
  if (!is.na(ref) && ref %in% levels(x)) {
    relevel(x, ref = ref)
  } else {
    x
  }
}

make_ward_3level <- function(data, ward_col = "Ward") {
  data %>%
    mutate(
      Ward_ED_ICU = case_when(
        .data[[ward_col]] == "Emergency" ~ "ED",
        .data[[ward_col]] == "ICU" ~ "ICU",
        TRUE ~ "Other"
      )
    )
}

# ============================================================
# 3) Phenotype stratum helpers
# ============================================================

classify_from_index <- function(data, index_cols, resistant_label, susceptible_label,
                                treat_i_as_r = FALSE) {
  idx_present <- intersect(index_cols, names(data))
  if (length(idx_present) == 0) {
    stop(
      "None of the index columns were found in the dataset. Looked for: ",
      paste(index_cols, collapse = ", ")
    )
  }

  bucket_matrix <- sapply(idx_present, function(col) ast_bucket(data[[col]]))
  if (is.null(dim(bucket_matrix))) {
    bucket_matrix <- matrix(bucket_matrix, ncol = 1)
    colnames(bucket_matrix) <- idx_present
  }

  any_r  <- apply(bucket_matrix, 1, function(v) any(v %in% "R", na.rm = TRUE))
  any_i  <- apply(bucket_matrix, 1, function(v) any(v %in% "I", na.rm = TRUE))
  any_si <- apply(bucket_matrix, 1, function(v) any(v %in% c("S", "I"), na.rm = TRUE))

  if (isTRUE(treat_i_as_r)) {
    any_r <- any_r | any_i
  }

  out <- rep(NA_character_, nrow(data))
  out[any_r] <- resistant_label
  out[!any_r & any_si] <- susceptible_label
  out
}

add_3gc_group <- function(data, index_cols = c("CTX", "CAZ"), treat_i_as_r = FALSE) {
  data %>%
    mutate(
      g3_group = classify_from_index(
        data = cur_data_all(),
        index_cols = index_cols,
        resistant_label = "3GCR",
        susceptible_label = "3GCS",
        treat_i_as_r = treat_i_as_r
      )
    )
}

add_carb_group <- function(data, index_cols = c("IPM", "ETP", "MEM", "MERO"),
                           treat_i_as_r = FALSE) {
  data %>%
    mutate(
      carb_group = classify_from_index(
        data = cur_data_all(),
        index_cols = index_cols,
        resistant_label = "CR",
        susceptible_label = "CS",
        treat_i_as_r = treat_i_as_r
      )
    )
}

# ============================================================
# 4) Model helpers
# ============================================================

tidy_glm_wald_or <- function(model, conf.level = 0.95, exponentiate = TRUE) {
  sm <- summary(model)$coefficients

  out <- tibble(
    term = rownames(sm),
    estimate = sm[, "Estimate"],
    std.error = sm[, "Std. Error"],
    statistic = sm[, "z value"],
    p.value = sm[, "Pr(>|z|)"]
  )

  z <- qnorm(1 - (1 - conf.level) / 2)
  out <- out %>%
    mutate(
      conf.low = estimate - z * std.error,
      conf.high = estimate + z * std.error
    )

  if (isTRUE(exponentiate)) {
    out <- out %>%
      mutate(
        estimate = exp(estimate),
        conf.low = exp(conf.low),
        conf.high = exp(conf.high)
      )
  }

  out
}

max_vif_safe <- function(model) {
  aliased_terms <- tryCatch(alias(model)$Complete, error = function(e) NULL)
  if (!is.null(aliased_terms) && any(aliased_terms)) {
    return(NA_real_)
  }

  tryCatch({
    vif_out <- car::vif(model, type = "predictor")
    max(as.numeric(vif_out), na.rm = TRUE)
  }, error = function(e) NA_real_)
}

fit_binary_resistance_model <- function(data, outcome_col,
                                        formula = Resistance ~ Year + Specimen_type + Hospital * Ward_ED_ICU,
                                        count_i_as_r_outcome = FALSE,
                                        conf.level = 0.95) {
  tmp <- data %>%
    mutate(bucket = ast_bucket(.data[[outcome_col]])) %>%
    filter(!is.na(bucket)) %>%
    mutate(
      Resistance = if (count_i_as_r_outcome) {
        as.integer(bucket %in% c("R", "I"))
      } else {
        as.integer(bucket == "R")
      }
    ) %>%
    drop_na(Resistance, Year, Specimen_type, Hospital, Ward_ED_ICU)

  if (nrow(tmp) == 0) {
    return(list(status = "skipped", reason = "no usable AST after filtering", data = tmp))
  }

  if (dplyr::n_distinct(tmp$Resistance) < 2L) {
    return(list(status = "skipped", reason = "outcome constant within stratum", data = tmp))
  }

  warnings_seen <- character(0)

  model <- tryCatch(
    withCallingHandlers(
      glm(formula = formula, family = binomial, data = tmp),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if (is.null(model)) {
    return(list(
      status = "skipped",
      reason = "glm failed to fit",
      data = tmp,
      warnings = unique(warnings_seen)
    ))
  }

  summary_tbl <- tidy_glm_wald_or(model, conf.level = conf.level, exponentiate = TRUE)

  global_tbl <- tryCatch({
    drop1(model, test = "Chisq") %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(Variable != "<none>") %>%
      transmute(Variable, Global_p = `Pr(>Chi)`)
  }, error = function(e) {
    tibble(Variable = NA_character_, Global_p = NA_real_)
  })

  metrics_tbl <- tibble(
    N = nrow(tmp),
    Events = sum(tmp$Resistance),
    AIC = AIC(model),
    Max_VIF = max_vif_safe(model),
    warnings_n = length(unique(warnings_seen))
  )

  list(
    status = "fit",
    reason = NA_character_,
    data = tmp,
    model = model,
    summary = summary_tbl,
    global = global_tbl,
    metrics = metrics_tbl,
    warnings = unique(warnings_seen)
  )
}

# ============================================================
# 5) Main model runner for GitHub use
# ============================================================

run_models_by_stratum <- function(data,
                                  abx_list,
                                  stratum = c("all", "3GCR", "3GCS", "CR", "CS"),
                                  index_cols = NULL,
                                  treat_i_as_r_for_index = FALSE,
                                  count_i_as_r_outcome = FALSE,
                                  ref_hospital = "HCN",
                                  ref_specimen = "Blood",
                                  ref_ward = "ICU",
                                  formula = Resistance ~ Year + Specimen_type + Hospital * Ward_ED_ICU,
                                  bug_name = NA_character_) {
  stratum <- match.arg(stratum)

  required_covars <- c("Hospital", "Specimen_type", "Ward", "Year")
  missing_covars <- setdiff(required_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Missing required covariate columns: ", paste(missing_covars, collapse = ", "))
  }

  dat <- data %>%
    make_ward_3level() %>%
    mutate(
      Hospital = safe_relevel(Hospital, ref_hospital),
      Specimen_type = safe_relevel(Specimen_type, ref_specimen),
      Ward_ED_ICU = safe_relevel(Ward_ED_ICU, ref_ward),
      Year = as.numeric(Year)
    )

  if (stratum %in% c("3GCR", "3GCS")) {
    if (is.null(index_cols)) index_cols <- c("CTX", "CAZ")
    dat <- add_3gc_group(dat, index_cols = index_cols, treat_i_as_r = treat_i_as_r_for_index) %>%
      filter(g3_group == stratum)
  }

  if (stratum %in% c("CR", "CS")) {
    if (is.null(index_cols)) index_cols <- c("IPM", "ETP", "MEM", "MERO")
    dat <- add_carb_group(dat, index_cols = index_cols, treat_i_as_r = treat_i_as_r_for_index) %>%
      filter(carb_group == stratum)
  }

  dat <- dat %>% drop_na(Year, Specimen_type, Hospital, Ward_ED_ICU)
  if (nrow(dat) == 0) {
    stop("After filtering to stratum '", stratum, "', there are 0 rows.")
  }

  abx_present <- intersect(abx_list, names(dat))
  if (length(abx_present) == 0) {
    stop("None of the requested antibiotic columns were found in the dataset.")
  }

  # Exclude index drugs when they define the stratum.
  if (stratum != "all" && !is.null(index_cols)) {
    abx_present <- setdiff(abx_present, intersect(abx_present, index_cols))
  }
  if (length(abx_present) == 0) {
    stop("No antibiotics left to model after excluding index columns.")
  }

  summary_list <- list()
  global_list <- list()
  metrics_list <- list()
  log_list <- list()
  model_list <- list()

  for (abx in abx_present) {
    fit <- fit_binary_resistance_model(
      data = dat,
      outcome_col = abx,
      formula = formula,
      count_i_as_r_outcome = count_i_as_r_outcome
    )

    log_list[[abx]] <- tibble(
      bug = bug_name,
      Stratum = stratum,
      Antibiotic = abx,
      status = fit$status,
      reason = fit$reason %||% NA_character_,
      N = nrow(fit$data),
      Events = if ("Resistance" %in% names(fit$data)) sum(fit$data$Resistance) else 0,
      warnings_n = length(fit$warnings %||% character(0))
    )

    if (fit$status != "fit") next

    summary_list[[abx]] <- fit$summary %>%
      mutate(
        bug = bug_name,
        Stratum = stratum,
        Antibiotic = abx,
        warnings_n = length(fit$warnings)
      )

    global_list[[abx]] <- fit$global %>%
      mutate(bug = bug_name, Stratum = stratum, Antibiotic = abx)

    metrics_list[[abx]] <- fit$metrics %>%
      mutate(bug = bug_name, Stratum = stratum, Antibiotic = abx, .before = 1)

    model_list[[abx]] <- fit$model
  }

  summary_df <- bind_rows(summary_list)
  global_df <- bind_rows(global_list)
  metrics_df <- bind_rows(metrics_list)
  log_df <- bind_rows(log_list)

  year_effect_df <- summary_df %>%
    filter(term == "Year") %>%
    transmute(
      bug,
      Stratum,
      Antibiotic,
      aOR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value,
      N = metrics_df$N[match(Antibiotic, metrics_df$Antibiotic)],
      Events = metrics_df$Events[match(Antibiotic, metrics_df$Antibiotic)]
    )

  list(
    summary = summary_df,
    year_effect = year_effect_df,
    global = global_df,
    metrics = metrics_df,
    log = log_df,
    models = model_list
  )
}

export_model_results <- function(results, prefix) {
  write_xlsx(results$summary,     paste0(prefix, "_summary.xlsx"))
  write_xlsx(results$year_effect, paste0(prefix, "_year_effect.xlsx"))
  write_xlsx(results$global,      paste0(prefix, "_global_pvalues.xlsx"))
  write_xlsx(results$metrics,     paste0(prefix, "_metrics.xlsx"))
  write_xlsx(results$log,         paste0(prefix, "_runlog.xlsx"))
}

# ============================================================
# 6) Plotting helper: yearly trends by species and carbapenem status
# ============================================================

build_yearly_trend_plot <- function(xlsx_path,
                                    sheet_yearly = "yearly_counts",
                                    sheet_trend = "trend_summary",
                                    p_col = "p_model",
                                    alpha = 0.05,
                                    drug_codes = c("amk", "cip", "ipm", "etp", "atm", "sxt", "gen", "tgc", "fof"),
                                    bug_order = c("K. pneumoniae", "E. coli", "Acinetobacter spp.", "P. aeruginosa")) {
  tol_muted9 <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499"
  )

  drug_cols <- setNames(tol_muted9[seq_along(drug_codes)], drug_codes)
  drug_labs <- setNames(drug_label(drug_codes), drug_codes)

  yearly <- read_excel(xlsx_path, sheet = sheet_yearly)
  trend  <- read_excel(xlsx_path, sheet = sheet_trend)

  if (!"prop" %in% names(yearly) && all(c("n_res", "n_tested") %in% names(yearly))) {
    yearly <- yearly %>% mutate(prop = n_res / n_tested)
  }

  plot_df <- yearly %>%
    mutate(
      year = as.integer(year),
      bug = toupper(as.character(bug)),
      bug_lab = bug_label(bug),
      carb_group = factor(carb_group, levels = c("CS", "CR")),
      drug = tolower(as.character(drug))
    ) %>%
    filter(drug %in% drug_codes) %>%
    left_join(
      trend %>%
        mutate(
          bug = toupper(as.character(bug)),
          carb_group = factor(carb_group, levels = c("CS", "CR")),
          drug = tolower(as.character(drug))
        ) %>%
        select(bug, carb_group, drug, p_model, p_adj_fdr),
      by = c("bug", "carb_group", "drug")
    ) %>%
    mutate(
      p_use = .data[[p_col]],
      sig = if_else(!is.na(p_use) & p_use < alpha, "Significant", "Not significant"),
      sig = factor(sig, levels = c("Significant", "Not significant"))
    ) %>%
    group_by(bug, bug_lab, carb_group, drug) %>%
    complete(year = seq(min(year, na.rm = TRUE), max(year, na.rm = TRUE), by = 1)) %>%
    fill(sig, p_use, p_model, p_adj_fdr, .direction = "downup") %>%
    ungroup() %>%
    mutate(
      n_tested = as.integer(n_tested),
      n_res = as.integer(n_res),
      prop = as.numeric(prop),
      prop_plot = if_else(!is.na(n_tested) & n_tested > 0, prop, NA_real_),
      bug_lab = factor(bug_lab, levels = bug_order),
      carb_group = factor(carb_group, levels = c("CS", "CR")),
      drug = factor(drug, levels = drug_codes),
      sig = factor(as.character(sig), levels = c("Significant", "Not significant"))
    ) %>%
    filter(!is.na(prop_plot)) %>%
    droplevels()

  year_breaks <- sort(unique(plot_df$year))

  ggplot(plot_df, aes(x = year, y = prop_plot)) +
    geom_area(
      aes(fill = drug, alpha = sig, group = drug),
      position = "identity",
      na.rm = TRUE
    ) +
    geom_line(
      aes(colour = drug, linetype = sig, group = drug),
      linewidth = 0.85,
      na.rm = TRUE
    ) +
    facet_grid(
      bug_lab ~ carb_group,
      labeller = labeller(
        carb_group = c(CS = "Carbapenem susceptible", CR = "Carbapenem resistant")
      )
    ) +
    scale_fill_manual(values = drug_cols, breaks = drug_codes, labels = drug_labs, drop = FALSE) +
    scale_colour_manual(values = drug_cols, breaks = drug_codes, labels = drug_labs, drop = FALSE, guide = "none") +
    scale_linetype_manual(values = c(Significant = "solid", `Not significant` = "dashed")) +
    scale_alpha_manual(values = c(Significant = 0.28, `Not significant` = 0.14), guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_x_continuous(breaks = year_breaks) +
    labs(x = "Year", y = "Resistance (%)", fill = NULL, linetype = "Model M2") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = list(alpha = 1)),
      linetype = guide_legend(order = 2)
    )
}

# ============================================================
# 7) Nested model comparison helper
# ============================================================

compare_nested_models <- function(data,
                                  abx_list,
                                  formula_small,
                                  formula_large,
                                  ref_hospital = "HCN",
                                  ref_specimen = "Blood",
                                  ref_ward = "ICU") {
  dat <- data %>%
    make_ward_3level() %>%
    mutate(
      Hospital = safe_relevel(Hospital, ref_hospital),
      Specimen_type = safe_relevel(Specimen_type, ref_specimen),
      Ward = safe_relevel(Ward, ref_ward),
      Ward_ED_ICU = safe_relevel(Ward_ED_ICU, ref_ward),
      Year = as.numeric(Year)
    )

  out <- list()

  for (abx in intersect(abx_list, names(dat))) {
    tmp <- dat %>%
      mutate(bucket = ast_bucket(.data[[abx]])) %>%
      filter(!is.na(bucket)) %>%
      mutate(Resistance = as.integer(bucket == "R")) %>%
      drop_na(Resistance, Year, Hospital, Specimen_type)

    if (dplyr::n_distinct(tmp$Resistance) < 2L) next

    model_small <- glm(formula_small, family = binomial, data = tmp)
    model_large <- glm(formula_large, family = binomial, data = tmp)
    lrt <- anova(model_small, model_large, test = "Chisq")

    out[[abx]] <- tibble(
      Antibiotic = abx,
      AIC_small = AIC(model_small),
      AIC_large = AIC(model_large),
      LRT_p_value = lrt$`Pr(>Chi)`[2]
    )
  }

  bind_rows(out)
}

# ============================================================
# 8) Example usage
# ============================================================

# Example antibiotic panels
abx_panels <- list(
  ECOL = c("CTX", "CAZ", "FEP", "GEN", "CIP", "AMP", "SXT", "IPM", "ETP", "AMK", "TGC", "ATM", "TZP"),
  KPN  = c("CTX", "CAZ", "FEP", "IPM", "ETP", "AMK", "GEN", "CIP", "ATM", "SXT", "TZP"),
  PSA  = c("CAZ", "FEP", "CIP", "AMK", "IPM", "TZP", "ATM", "MERO"),
  ACN  = c("SXT", "IPM", "AMK", "CIP", "GEN")
)

carb_index_by_bug <- list(
  ECOL = c("IPM", "ETP", "MEM", "MERO"),
  KPN  = c("IPM", "ETP", "MEM", "MERO"),
  PSA  = c("IPM", "MEM", "MERO"),
  ACN  = c("IPM", "MEM", "MERO")
)

# Example runs (uncomment after loading your objects, e.g. ECOL / KPN / PSA / ACN)
# res_kpn_3gcr <- run_models_by_stratum(
#   data = KPN,
#   abx_list = abx_panels$KPN,
#   stratum = "3GCR",
#   index_cols = c("CTX", "CAZ"),
#   bug_name = "KPN"
# )
# export_model_results(res_kpn_3gcr, "KPN_M2_within3GCR")
#
# res_ecol_cr <- run_models_by_stratum(
#   data = ECOL,
#   abx_list = abx_panels$ECOL,
#   stratum = "CR",
#   index_cols = carb_index_by_bug$ECOL,
#   bug_name = "ECOL"
# )
# export_model_results(res_ecol_cr, "ECOL_M2_withinCR")
#
# p <- build_yearly_trend_plot("trend_test_CR_vs_CS_by_bug_drug_ANALYSIS_ONLY.xlsx")
# ggsave("JAC_area_CS_vs_CR_yearly_trends.pdf", p, width = 12, height = 9)
# ggsave("JAC_area_CS_vs_CR_yearly_trends.png", p, width = 12, height = 9, dpi = 300)
