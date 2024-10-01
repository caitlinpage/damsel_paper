library(Damsel)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
damsel_counts <- readRDS("data/damsel_counts.rds")
gatc_regions <- getGatcRegions(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)$regions

damsel_dge <- makeDGE(damsel_counts, min.samples = 2)

damsel_dm <- testDmRegions(damsel_dge, gatc_regions)

damsel_peaks <- identifyPeaks(damsel_dm)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
genes <- collateGenes(genes = txdb, regions = gatc_regions, org.Db = org.Dm.eg.db)

damsel_genes <- annotatePeaksGenes(damsel_peaks, genes, gatc_regions)

saveRDS(damsel_dm, "output/damsel_dm.rds")
saveRDS(damsel_peaks, "output/damsel_peaks.rds")
saveRDS(damsel_genes, "output/damsel_genes.rds")
