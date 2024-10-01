
library(Damsel)
library(plyranges)
library(dplyr)
library(tidyr)
library(edgeR)

library(BSgenome.Dmelanogaster.UCSC.dm6)

gatc_regions <- getGatcRegions(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)$regions
damsel_counts <- readRDS("data/damsel_counts.rds")


matrix <- as.matrix(damsel_counts[, grepl("bam", colnames(damsel_counts), ignore.case = TRUE)])
rownames(matrix) <- damsel_counts$Position
group = c("Dam", "Sd", "Dam", "Sd")
design = model.matrix(~group)
y = DGEList(matrix, group = group)
keep <- rowSums(cpm(y)>=0.5) >= 2
y = y[keep, ,keep.lib.sizes=FALSE]
y = calcNormFactors(y)
y = estimateDisp(y, robust = T, design = design)

fit = glmFit(y, design = design)
lrt = glmLRT(fit, coef=2)
de.Sd <- decideTestsDGE(lrt, lfc = 1)
lrt$table$significant <- de.Sd

vissers_dm <- data.frame(lrt$table)
vissers_dm$significant <- data.frame(de.Sd)$groupSd

write.table(vissers_dm, file='output/lrt_sd.txt', quote=F)
write.table(keep, file='output/keep', quote=F, col.names = FALSE)

system2("python3", args=c("code/call_peaks.py",
        "output/keep", "output/lrt_sd.txt", ">",
        "output/vissers_peaks.txt"))
vissers_peaks <- read.table("output/vissers_peaks.txt")

names(vissers_peaks) <- c('seqnames', 'start', 'end', "tags", 'pen', 'aveLogFC', 'sig')
vissers_peaks <- vissers_peaks %>% mutate(big = ifelse(tags > 2 | pen == 1, TRUE, FALSE))

vissers_peaks_mod <- vissers_peaks
names(vissers_peaks_mod) <- c('seqnames', 'start', 'end', "tags", 'pen', 'aveLogFC', 'sig', "big")


peaks_per_chr <- vissers_peaks_mod %>%
  mutate(end_chr = ifelse(start > lead(start), TRUE, NA),
         row_num = 1:nrow(.),
         end_chr = ifelse(row_num == nrow(.), TRUE, end_chr),
         end_pos = ifelse(end_chr == TRUE, row_num, NA)) %>%
  fill(end_pos, .direction = "up") %>%
  distinct(end_pos) %>%
  mutate(num_per_chr = end_pos - lag(end_pos),
         num_per_chr = ifelse(is.na(num_per_chr), end_pos, num_per_chr)) %>%
  .$num_per_chr

chromosomes <- sub("-.*", "", rownames(vissers_dm)) %>% unique()

if(length(chromosomes) == length(peaks_per_chr)) {

  chr_vec <- character()
  for(i in 1:length(peaks_per_chr)) {
    chr_vec <- c(chr_vec, rep(chromosomes[i], times = peaks_per_chr[i]))
  }

  vissers_peaks_mod$seqnames <- chr_vec
}

saveRDS(vissers_dm, "output/vissers_dm.rds")
saveRDS(vissers_peaks, "output/vissers_peaks.rds")
saveRDS(vissers_peaks_mod, "output/vissers_peaks_mod.rds")

