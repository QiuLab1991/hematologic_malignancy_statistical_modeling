suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
})

set.seed(2025)

exposure_dir <- "data/pQTL_exposures"
outcome_dir <- "data/gwas_outcomes"
output_dir <- "results/genetic_mr"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

exposure_files <- list.files(
  exposure_dir,
  pattern = "\\.csv(\\.gz)?$",
  full.names = TRUE
)

outcome_files <- list.files(
  outcome_dir,
  pattern = "\\.csv(\\.gz)?$",
  full.names = TRUE
)

if (length(exposure_files) == 0) stop("No pQTL exposure files were found.")
if (length(outcome_files) == 0) stop("No GWAS outcome files were found.")

all_mr <- list()
all_heterogeneity <- list()
all_pleiotropy <- list()
all_harmonised <- list()

for (exposure_file in exposure_files) {
  source_id <- sub("\\.csv(\\.gz)?$", "", basename(exposure_file))

  exposure_dat <- read_exposure_data(
    filename = exposure_file,
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    pval_col = "pval",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf.exposure",
    phenotype_col = "gene",
    id_col = "id.exposure",
    samplesize_col = "samplesize.exposure",
    chr_col = "CHR",
    pos_col = "BP",
    clump = FALSE
  )

  for (outcome_file in outcome_files) {
    outcome_id <- sub("\\.csv(\\.gz)?$", "", basename(outcome_file))
    outcome_raw <- fread(outcome_file)

    required_outcome_columns <- c("SNP", "beta", "se", "P", "A1", "A2", "MAF")
    missing_outcome_columns <- setdiff(required_outcome_columns, names(outcome_raw))
    if (length(missing_outcome_columns) > 0) {
      stop("Missing outcome columns in ", basename(outcome_file), ": ",
           paste(missing_outcome_columns, collapse = ", "))
    }

    if (!"samplesize" %in% names(outcome_raw)) {
      if (all(c("nCase", "nControl") %in% names(outcome_raw))) {
        outcome_raw[, samplesize := nCase + nControl]
      } else {
        stop("Outcome file must contain either samplesize or both nCase and nControl: ",
             basename(outcome_file))
      }
    }

    outcome_raw <- outcome_raw[SNP %in% unique(exposure_dat$SNP)]
    if (nrow(outcome_raw) == 0) next

    outcome_raw[, outcome := outcome_id]
    outcome_raw[, id := outcome_id]

    outcome_dat <- format_data(
      dat = as.data.frame(outcome_raw),
      type = "outcome",
      phenotype_col = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "P",
      samplesize_col = "samplesize",
      id_col = "id"
    )

    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    dat <- dat[dat$mr_keep, , drop = FALSE]

    if (nrow(dat) == 0) next

    dat$Source <- source_id
    key <- paste(source_id, outcome_id, sep = "__")
    all_harmonised[[key]] <- dat

    mr_res <- mr(
      dat,
      method_list = c("mr_wald_ratio", "mr_ivw")
    )

    if (nrow(mr_res) > 0) {
      mr_res <- generate_odds_ratios(mr_res)
      mr_res$Source <- source_id
      all_mr[[key]] <- mr_res
    }

    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    if (!is.null(het) && nrow(het) > 0) {
      het$Source <- source_id
      all_heterogeneity[[key]] <- het
    }

    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    if (!is.null(pleio) && nrow(pleio) > 0) {
      pleio$Source <- source_id
      all_pleiotropy[[key]] <- pleio
    }
  }
}

mr_results <- bind_rows(all_mr)
if (nrow(mr_results) == 0) stop("No MR results were generated.")

mr_results <- mr_results %>%
  group_by(Source, id.outcome) %>%
  mutate(FDR = p.adjust(pval, method = "BH")) %>%
  ungroup()

fwrite(as.data.table(mr_results), file.path(output_dir, "mr_results.csv"))

if (length(all_harmonised) > 0) {
  harmonised <- bind_rows(all_harmonised)
  fwrite(as.data.table(harmonised), file.path(output_dir, "harmonised_instruments.csv"))
}

if (length(all_heterogeneity) > 0) {
  heterogeneity <- bind_rows(all_heterogeneity)
  fwrite(as.data.table(heterogeneity), file.path(output_dir, "mr_heterogeneity.csv"))
}

if (length(all_pleiotropy) > 0) {
  pleiotropy <- bind_rows(all_pleiotropy)
  fwrite(as.data.table(pleiotropy), file.path(output_dir, "mr_egger_intercept.csv"))
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))




