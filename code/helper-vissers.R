
library(Damsel)
library(plyranges)
library(dplyr)
library(edgeR)

test_vissers_dm <- function(damsel_counts, modify_group=NULL) {
  matrix <- as.matrix(damsel_counts[, grepl("bam", colnames(damsel_counts), ignore.case = TRUE)])
  rownames(matrix) <- damsel_counts$Position
  #my code
  if(is.null(modify_group)) {
    group <- c("Dam", "Sd", "Dam", "Sd")
  } else {
    group <- modify_group
  }
  design = model.matrix(~group)
  y = DGEList(matrix, group = group)
  keep <- rowSums(cpm(y)>=0.5) >= 2
  y = y[keep, ,keep.lib.sizes=FALSE]
  y = calcNormFactors(y)
  y = estimateDisp(y, robust = T, design = design)

  fit = glmFit(y, design = design)
  lrt = glmLRT(fit, coef=2)
  #code uses decideTestsDGE - now deprecated
  #replaced with decideTests
  de.Sd <- decideTests(lrt, lfc = 1)
  #lrt$table$significant <- de.Sd

  vissers_dm <- data.frame(lrt$table)
  vissers_dm$significant <- data.frame(de.Sd)[,1]
  list(design=design, results=vissers_dm)
}


test_vissers_peaks <- function(damsel_counts, modify_group=NULL) {
  matrix <- as.matrix(damsel_counts[, grepl("bam", colnames(damsel_counts), ignore.case = TRUE)])
  rownames(matrix) <- damsel_counts$Position
  #my code
  if(is.null(modify_group)) {
    group <- c("Dam", "Sd", "Dam", "Sd")
  } else {
    group <- modify_group
  }
  design = model.matrix(~group)
  y = DGEList(matrix, group = group)
  keep <- rowSums(cpm(y)>=0.5) >= 2
  y = y[keep, ,keep.lib.sizes=FALSE]
  y = calcNormFactors(y)
  y = estimateDisp(y, robust = T, design = design)

  fit = glmFit(y, design = design)
  lrt = glmLRT(fit, coef=2)
  de.Sd <- decideTests(lrt, lfc = 1)
  lrt$table$significant <- de.Sd

  vissers_dm <- data.frame(lrt$table)
  vissers_dm$significant <- data.frame(de.Sd)[,1]

  write.table(vissers_dm, file='../output/lrt_sd.txt', quote=F)
  write.table(keep, file='../output/keep', quote=F, col.names = FALSE)

  system2("python3", args=c("../code/call_peaks.py",
                            "../output/keep", "../output/lrt_sd.txt", ">",
                            "../output/fp_vissers_peaks.txt"))
  vissers_peaks <- read.table("../output/fp_vissers_peaks.txt")

  names(vissers_peaks) <- c('seqnames', 'start', 'end', "tags", 'pen', 'aveLogFC', 'sig')
  vissers_peaks
}
