suppressPackageStartupMessages({
  library(data.table)
  library(Biobase)
  library(Mfuzz)
})

set.seed(2025)

input_file <- "data/prot.dia.Blood.csv.gz"
candidate_file <- "data/cox_prots.txt"
output_dir <- "results/pseudotime"
event_var <- "Blood"
time_var <- "HR.Blood"
n_clusters <- 4
maximum_years_to_diagnosis <- 15
bin_width_years <- 3

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)
candidate_table <- fread(candidate_file, header = TRUE)
if (ncol(candidate_table) < 1) stop("candidate_file must contain at least one column.")
candidate_proteins <- unique(toupper(as.character(candidate_table[[1]])))

required_columns <- c(event_var, time_var)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

name_map <- setNames(names(dat), toupper(names(dat)))
matched_proteins <- unname(name_map[intersect(candidate_proteins, names(name_map))])
matched_proteins <- matched_proteins[!is.na(matched_proteins)]
if (length(matched_proteins) < n_clusters) stop("Too few candidate proteins were found in the input data.")

cases <- dat[
  get(event_var) == 1 &
    !is.na(get(time_var)) &
    get(time_var) > 0 &
    get(time_var) < maximum_years_to_diagnosis
]

if (nrow(cases) == 0) stop("No incident cases were available in the requested time window.")

expression_matrix <- as.matrix(cases[, ..matched_proteins])
storage.mode(expression_matrix) <- "double"

keep_nonempty <- colSums(!is.na(expression_matrix)) > 0
expression_matrix <- expression_matrix[, keep_nonempty, drop = FALSE]

for (j in seq_len(ncol(expression_matrix))) {
  x <- expression_matrix[, j]
  med <- median(x, na.rm = TRUE)
  x[is.na(x)] <- med
  expression_matrix[, j] <- x
}

keep_variable <- apply(expression_matrix, 2, sd) > 0
expression_matrix <- expression_matrix[, keep_variable, drop = FALSE]

pseudo_time <- -cases[[time_var]]
bin_breaks <- seq(-maximum_years_to_diagnosis, 0, by = bin_width_years)
if (tail(bin_breaks, 1) != 0) bin_breaks <- c(bin_breaks, 0)
bin_centers <- head(bin_breaks, -1) + diff(bin_breaks) / 2

time_bin <- cut(
  pseudo_time,
  breaks = bin_breaks,
  labels = bin_centers,
  include.lowest = TRUE,
  right = FALSE
)

grouped <- aggregate(
  expression_matrix,
  by = list(PseudoTime = as.numeric(as.character(time_bin))),
  FUN = mean
)
grouped <- grouped[order(grouped$PseudoTime), ]
rownames(grouped) <- grouped$PseudoTime
grouped$PseudoTime <- NULL

expr_t <- t(as.matrix(grouped))
if (nrow(expr_t) < n_clusters) stop("Too few proteins remained for the requested number of clusters.")
if (ncol(expr_t) < 2) stop("At least two pseudo-time bins are required.")

eset <- ExpressionSet(assayData = expr_t)
eset <- standardise(eset)
m <- mestimate(eset)
cl <- mfuzz(eset, centers = n_clusters, m = m)

membership <- as.data.table(cl$membership, keep.rownames = "Protein")
setnames(membership, 2:ncol(membership), paste0("Cluster", seq_len(n_clusters), "_membership"))
cluster_assignment <- data.table(
  Protein = names(cl$cluster),
  AssignedCluster = as.integer(cl$cluster)
)
membership <- merge(membership, cluster_assignment, by = "Protein", all.x = TRUE, sort = FALSE)
membership[, MaxMembership := apply(.SD, 1, max), .SDcols = patterns("_membership$")]

group_means <- as.data.table(expr_t, keep.rownames = "Protein")

fwrite(membership, file.path(output_dir, "mfuzz_membership.csv"))
fwrite(group_means, file.path(output_dir, "pseudo_time_group_means.csv"))

pdf(file.path(output_dir, "mfuzz_clusters.pdf"), width = 6, height = 6)
mfuzz.plot(
  eset,
  cl = cl,
  mfrow = c(2, 2),
  new.window = FALSE,
  time.labels = colnames(expr_t)
)
dev.off()

saveRDS(list(eset = eset, clustering = cl, fuzzification = m), file.path(output_dir, "mfuzz_model.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
