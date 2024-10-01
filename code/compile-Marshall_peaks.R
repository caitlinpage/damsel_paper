library(rtracklayer)
library(plyranges)
library(dplyr)
library(tidyr)
marshall_peaks_1 <- rtracklayer::import("data/sd_1_SRR794884-vs-Dam.kde-norm.gatc-FDR0.01.peaks.gff")
marshall_peaks_2 <- rtracklayer::import("data/sd_2_SRR7948877-vs-Dam.kde-norm.gatc-FDR0.01.peaks.gff")



marshall_union <- union_ranges(marshall_peaks_1, marshall_peaks_2) %>% data.frame() %>%
  mutate(num = 1:n()) %>% as_granges()

marshall_overlap1 <- pair_overlaps(marshall_union, marshall_peaks_1) %>% data.frame()
marshall_overlap2 <- pair_overlaps(marshall_union, marshall_peaks_2) %>% data.frame()

marshall_union <- marshall_union %>%
  data.frame() %>%
  mutate(rep1 = marshall_overlap1[match(.$num, marshall_overlap1$num), "FDR"],
         rep2 = marshall_overlap2[match(.$num, marshall_overlap2$num), "FDR"],
         rep1 = as.double(rep1), rep2 = as.double(rep2)) %>%
  mutate(FDR = case_when(is.na(rep1) ~ rep2, is.na(rep2) ~ rep1, rep1 <= rep2 ~ rep1, TRUE ~ rep2)) %>%
  replace(is.na(.), 1) %>%
  mutate(source = ifelse(FDR == rep1, "rep1", "rep2")) %>%
  mutate(score = ifelse(source == "rep1", marshall_overlap1[match(.$num, marshall_overlap1$num), "score"],
                        marshall_overlap2[match(.$num, marshall_overlap2$num), "score"])) %>%
  mutate(seqnames = paste0("chr", seqnames)) %>% .[,c(1:6,9,11)] %>% .[order(.$FDR),] %>% mutate(rank_fdr = 1:n())
marshall_union

saveRDS(marshall_union, "output/marshall_peaks.rds")
