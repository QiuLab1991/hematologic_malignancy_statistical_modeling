suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

set.seed(2025)

input_file <- "data/prot.dia.Blood.csv.gz"
output_dir <- "results/proteomics"
event_var <- "Blood"
time_var <- "HR.Blood"
cancer_var <- "Allcancer"
region_var <- "Assessment_distinct"
protein_columns <- 16:2935

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)

required_base <- c(event_var, time_var, cancer_var, region_var)
missing_base <- setdiff(required_base, names(dat))
if (length(missing_base) > 0) stop("Missing required columns: ", paste(missing_base, collapse = ", "))

dat <- dat[!is.na(get(time_var))]
dat <- dat[get(event_var) == 1 | (get(cancer_var) == 0 & get(event_var) == 0)]

drop_proteins <- intersect(c("GLIPR1", "NPM1", "PCOLCE"), names(dat))
if (length(drop_proteins) > 0) dat[, (drop_proteins) := NULL]

if (max(protein_columns) > ncol(dat)) stop("protein_columns exceeds the number of columns in the input data.")
protein_names <- names(dat)[protein_columns]
protein_names <- protein_names[vapply(dat[, ..protein_names], is.numeric, logical(1))]
if (length(protein_names) == 0) stop("No numeric protein columns were identified.")

categorical_vars <- intersect(
  c("Sex", "Ethnic_white", "Assessment_centre", "Assessment_distinct", "education",
    "Employment", "Household_income", "TDI", "BMI", "Smoking", "Drinking",
    "physical_activity", "Sleep_during", "HBP", "DM", "HPL", "CVD", "CKD"),
  names(dat)
)
dat[, (categorical_vars) := lapply(.SD, factor), .SDcols = categorical_vars]

recruitment_season_var <- if ("Recruitment_season" %in% names(dat)) {
  "Recruitment_season"
} else if ("Recuit_season" %in% names(dat)) {
  "Recuit_season"
} else {
  stop("Recruitment season column was not found.")
}

model_covariates <- list(
  M0 = character(0),
  M1 = c("Age", "Sex", "Ethnic_white", "BMI"),
  M2 = c("Age", "Sex", "Ethnic_white", "BMI", "education", "Household_income",
         "TDI", "physical_activity", "sample_age", recruitment_season_var),
  M3 = c("Age", "Sex", "Ethnic_white", "BMI", "education", "Household_income",
         "TDI", "physical_activity", "sample_age", recruitment_season_var,
         "HBP", "DM", "HPL", "CVD", "CKD"),
  M4 = c("Age", "Sex", "Ethnic_white", "BMI", "education", "Household_income",
         "TDI", "physical_activity", "sample_age", recruitment_season_var,
         "HBP", "DM", "HPL", "CVD", "CKD", "Assessment_centre", "Employment",
         "Smoking", "Drinking", "Diet_Quality_Score", "Sleep_during")
)

missing_covariates <- setdiff(unique(unlist(model_covariates)), names(dat))
if (length(missing_covariates) > 0) stop("Missing covariates: ", paste(missing_covariates, collapse = ", "))

fit_one_protein <- function(data, protein, covariates) {
  model_data <- data[, c(time_var, event_var, covariates), with = FALSE]
  model_data[, Protein := data[[protein]]]
  rhs <- c("Protein", covariates)
  form <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ", paste(rhs, collapse = " + "))
  )
  fit <- coxph(form, data = model_data, ties = "efron")
  sm <- summary(fit)
  beta <- unname(coef(fit)["Protein"])
  se <- unname(sm$coefficients["Protein", "se(coef)"])
  p <- unname(sm$coefficients["Protein", "Pr(>|z|)"])
  data.table(
    Protein = protein,
    Beta = beta,
    SE = se,
    HR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    P = p,
    N = fit$n,
    Events = fit$nevent
  )
}

run_screen <- function(data, label, models) {
  out <- rbindlist(lapply(models, function(model_name) {
    covariates <- model_covariates[[model_name]]
    res <- rbindlist(lapply(protein_names, function(protein) {
      tryCatch(
        fit_one_protein(data, protein, covariates),
        error = function(e) data.table(
          Protein = protein, Beta = NA_real_, SE = NA_real_, HR = NA_real_,
          CI_lower = NA_real_, CI_upper = NA_real_, P = NA_real_,
          N = NA_integer_, Events = NA_integer_
        )
      )
    }), fill = TRUE)
    res[, Model := model_name]
    res[, FDR := p.adjust(P, method = "BH")]
    res
  }), fill = TRUE)
  out[, AnalysisSet := label]
  out[]
}

england <- dat[get(region_var) == "England"]
scotland_wales <- dat[get(region_var) %in% c("Scotland", "Wales")]

n_random60 <- floor(0.60 * nrow(dat))
random60_index <- sample.int(nrow(dat), size = n_random60, replace = FALSE)
random60 <- dat[random60_index]
random40 <- dat[-random60_index]

results <- rbindlist(list(
  run_screen(dat, "Full", c("M0", "M1", "M2", "M3", "M4")),
  run_screen(england, "Discovery_England", c("M0", "M1", "M2")),
  run_screen(scotland_wales, "InternalReplication_ScotlandWales", c("M0", "M1", "M2")),
  run_screen(random60, "Robustness_Random60", c("M0", "M1", "M2")),
  run_screen(random40, "Robustness_Random40", c("M0", "M1", "M2"))
), fill = TRUE)

setcolorder(results, c(
  "Protein", "AnalysisSet", "Model", "N", "Events",
  "Beta", "SE", "HR", "CI_lower", "CI_upper", "P", "FDR"
))

fwrite(results, file.path(output_dir, "prospective_proteomic_cox_results.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
