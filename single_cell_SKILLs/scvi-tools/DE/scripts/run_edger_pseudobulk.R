#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_edger_pseudobulk.R counts.csv.gz samples.csv results.csv")

suppressPackageStartupMessages(library(edgeR))
counts_file <- args[[1]]
samples_file <- args[[2]]
results_file <- args[[3]]

counts_df <- read.csv(counts_file, check.names = FALSE, row.names = 1)
samples <- read.csv(samples_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(colnames(counts_df), samples$sample_id)) stop("Counts and sample metadata are not aligned")

samples$condition <- factor(samples$condition, levels = c("control", "disease"))
design <- model.matrix(~ condition, data = samples)
y <- DGEList(counts = as.matrix(counts_df), samples = samples, group = samples$condition)
keep <- filterByExpr(y, design = design)
if (sum(keep) < 2) stop("Fewer than two genes passed filterByExpr")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "conditiondisease")
tab <- topTags(test, n = Inf, sort.by = "none")$table

cpm_values <- cpm(y, log = FALSE)
disease <- samples$condition == "disease"
tab$mean_cpm_disease <- rowMeans(cpm_values[, disease, drop = FALSE])
tab$mean_cpm_control <- rowMeans(cpm_values[, !disease, drop = FALSE])
tab$tested <- TRUE
tab$gene <- rownames(tab)
tab <- tab[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR",
               "mean_cpm_disease", "mean_cpm_control", "tested")]
write.csv(tab, results_file, row.names = FALSE)

