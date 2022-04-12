wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(DESeq2)
library(tximport)
library(tximeta)
library(AnnotationDbi)
library(data.table)
library(EnhancedVolcano)
library(msigdbr)
library(ComplexHeatmap)
library(DOSE)
library(enrichplot)
library(fgsea)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#######################
fp <- file.path(list.files("SALMON_QUANT/", full.names = T), "quant.sf")
fn <- list.files("SALMON_QUANT/") %>% gsub(pattern = "_quant", replacement = "")


colData <- data.frame(
  files = fp,
  names = fn
)


df.meta <- read_csv("GSE153348_parsed_meta_WithClinicalAnnotation.csv")


df.meta$names  <- paste0("SRR12097", seq(378, 428))

df.meta$achieved_cr <- df.meta$Best_Response %in% c("CR", "CRp") %>% factor(levels = c(TRUE, FALSE))

df.meta$timepoint <- gsub(df.meta$timepoint, pattern = "timepoint: ", replacement = "") %>% gsub(pattern = "\\??", replacement = "")

df.meta$id <- df.meta$Sample_title %>% 
  gsub(pattern = "RNA-seq_", replacement = "") %>% 
  gsub(pattern = "_BL", replacement = "") %>% 
  gsub(pattern = "_REL", replacement = "") %>%
  gsub(pattern = "\\??", replacement = "")

id.paired <- df.meta %>% group_by(id) %>% summarise(n = n()) %>% subset(n == 2) %>% .$id

df.meta$is_responder <- factor((df.meta$Best_Response %in% c("CR", "CRp", "MLFS", "PR", "HI")), levels = c("FALSE", "TRUE"))
##########################
##
colData <- left_join(colData, df.meta)

#########################
txi <- tximeta(coldata = colData)
gse <- summarizeToGene(txi)
dds <- DESeqDataSet(gse, design = ~achieved_cr + timepoint)
dds <- DESeq(dds)
#########################

counts <- counts(dds, normalized = TRUE) 
rownames(counts) <- rownames(counts) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
# counts %>% 
#   as.data.frame() %>% 
#   mutate(ENSEMBL = rownames(.)) %>% 
#   left_join(., map) %>% 
#   data.table() -> counts


#Take paired samples
colData_paired <- colData[colData$id %in% id.paired, ]
txi.paired <- tximeta(coldata = colData_paired)
gse.paired <- summarizeToGene(txi.paired)
dds.paired <- DESeqDataSet(gse.paired, design = ~achieved_cr + timepoint)
dds.paired <- DESeq(dds.paired)


#Take bl samples
colData_bl <- colData %>% subset(timepoint == "BL")
txi.bl <- tximeta(coldata = colData_bl)
gse.bl <- summarizeToGene(txi.bl)
dds.bl <- DESeqDataSet(gse.bl, design = ~is_responder)
dds.bl <- DESeq(dds.bl)



# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")


res.rel <- results(dds, contrast = c("timepoint", "REL", "BL"))
res.cr <- results(dds, contrast = c("achieved_cr", "FALSE", "TRUE"))

# res.rel <- lfcShrink(dds, coef = "timepoint_REL_vs_BL")
# res.cr <- lfcShrink(dds, coef = "achieved_cr_FALSE_vs_TRUE")
# 
# 
res.paired.rel <- results(dds.paired, contrast = c("timepoint", "REL", "BL"))


rownames(res.cr) <- rownames(res.cr) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
rownames(res.rel) <- rownames(res.rel) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
rownames(res.paired.rel ) <- rownames(res.paired.rel) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()


AnnotationDbi::select(org.Hs.eg.db,
                      keys = rownames(res.cr) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist(),
                      keytype = "ENSEMBL", 
                      columns = "SYMBOL") -> map

res.cr %>% 
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  left_join(., map) %>% 
  data.table() -> df.res.cr

res.rel %>% 
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  left_join(., map) %>% 
  data.table() -> df.res.rel


res.paired.rel %>% 
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  left_join(., map) %>% 
  data.table() -> df.res.paired.rel

# gene.stats.cr <- df.res.cr$stat
# names(gene.stats.cr) <- df.res.cr$ENTREZID
# gene.stats.cr <- gene.stats.cr[!is.na(gene.stats.cr)]
# 
# gene.stats.rel <- df.res.rel$stat
# names(gene.stats.rel) <- df.res.rel$ENTREZID
# gene.stats.rel <- gene.stats.rel[!is.na(gene.stats.rel)]

gene.stats.paired.rel <- df.res.paired.rel$stat
names(gene.stats.paired.rel) <- df.res.paired.rel$ENTREZID
gene.stats.paired.rel <- gene.stats.paired.rel[!is.na(gene.stats.paired.rel)]

gene.stats.cr <- df.res.cr$log2FoldChange
names(gene.stats.cr) <- df.res.cr$ENTREZID
gene.stats.cr <- gene.stats.cr[!is.na(gene.stats.cr)]

gene.stats.rel <- df.res.rel$log2FoldChange
names(gene.stats.rel) <- df.res.rel$ENTREZID
gene.stats.rel <- gene.stats.rel[!is.na(gene.stats.rel)]

####
hm <- msigdbr(category = "H")
c1 <- msigdbr(category = "C1")
c2 <- msigdbr(category = "C2")
c3 <- msigdbr(category = "C3")
c4 <- msigdbr(category = "C4")
c5 <- msigdbr(category = "C5")
c6 <- msigdbr(category = "C6")
c7 <- msigdbr(category = "C7")
c8 <- msigdbr(category = "C8")

gs.ls <- list(hm, c2, c3)
gs.ls %>% 
  lapply(function(x){
    indx <- grep(x$gs_name, pattern = "stat5", ignore.case = T)
    x[indx, ]
  }) -> gs.ls

gs <- do.call(rbind, gs.ls)
gs$gs_name <- factor(gs$gs_name)
gs.name <- levels(gs$gs_name)
gs %>% 
  group_by(gs_name) %>% 
  group_map(function(x, y){
    x$gene_symbol
  }) -> gs
names(gs) <- gs.name

rm(gs.ls)



gs.hox.ls <- list(hm, c1, c2, c3, c4, c5, c6, c7, c8)
gs.hox.ls %>% 
  lapply(function(x){
    indx <- grep(x$gs_name, pattern = "HOX", ignore.case = F)
    x[indx, ]
  }) -> gs.hox.ls
gs.hox <- do.call(rbind, gs.hox.ls)
gs.hox$gs_name <- factor(gs.hox$gs_name)
gs.hox.name <- levels(gs.hox$gs_name)
gs.hox %>% 
  group_by(gs_name) %>% 
  group_map(function(x, y){
    x$gene_symbol
  }) -> gs.hox
names(gs.hox) <- gs.hox.name

####
# gsea.cr <- fgseaSimple(gs, gene.stats.cr, 1E5)
# gsea.cr <- gsea.cr[order(gsea.cr$padj), ]
# 
# gsea.rel <- fgseaSimple(gs, gene.stats.rel, 1E4)
# gsea.rel <- gsea.rel[order(gsea.rel$padj), ]
# 
# gsea.paired.rel <- fgseaSimple(gs, gene.stats.paired.rel, 1E4)
# gsea.paired.rel <- gsea.paired.rel[order(gsea.paired.rel$padj), ]
# df.norm.count <- counts(dds, normalized = TRUE)
# df.norm.count %>% 
#   as.data.frame() %>% 
#   mutate(ENSEMBL = rownames(.) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
#   pivot_longer(cols = 1:(ncol(.)-1), names_to = "names", values_to = "norm.count") %>% 
#   left_join(., df.meta) %>% 
#   left_join(., map) -> df.norm.count.long
#write_csv(file = "normalized_count_LONG.csv", x = df.norm.count.long)




# PCA ---------------------------------------------------------------------

#dds.vst <- DESeq2::vst(dds, blind = FALSE)
#DESeq2::plotPCA(dds.vst, intgroup = "achieved_cr", ntop = 100)


##1. Do gene expression profiles at BL&REL form distinct clusters for pateitns who achieved CR?
colData_cr.true <- colData %>% subset(achieved_cr == TRUE) 
txi_cr.true <- tximeta(coldata = colData_cr.true)
gse_cr.true <- summarizeToGene(txi_cr.true)
dds_cr.true <- DESeqDataSet(gse_cr.true, design = ~timepoint)
dds_cr.true <- DESeq(dds_cr.true)
# vst_cr.true <- vst(dds_cr.true)
# DESeq2::plotPCA(vst_cr.true, intgroup = "timepoint", ntop = 100) 

##1-2. How about only IDH1 patients?
# colData_cr.true_midh1 <- colData %>% subset(achieved_cr == TRUE) %>% subset(IDH_mutation == "IDH1")
# txi_cr.true_midh1 <- tximeta(coldata = colData_cr.true_midh1)
# gse_cr.true_midh1 <- summarizeToGene(txi_cr.true_midh1)
# dds_cr.true_midh1 <- DESeqDataSet(gse_cr.true_midh1, design = ~timepoint)
# dds_cr.true_midh1 <- DESeq(dds_cr.true_midh1)
# vst_cr.true_midh1 <- vst(dds_cr.true_midh1)
# DESeq2::plotPCA(vst_cr.true_midh1, intgroup = "timepoint", ntop = 100) 

##2. Do gene expression profiles at baseline form distinct clusters based on clinical response?
colData_tp.bl <- colData %>% subset(timepoint == "BL") 
txi_tp.bl<- tximeta(coldata = colData_tp.bl)
gse_tp.bl<- summarizeToGene(txi_tp.bl)
dds_tp.bl<- DESeqDataSet(gse_tp.bl, design = ~achieved_cr + is_responder)
dds_tp.bl<- DESeq(dds_tp.bl)
#vst_tp.bl<- vst(dds_tp.bl)
# DESeq2::plotPCA(vst_tp.bl, intgroup = "Best_Response", ntop = 100) +  ###save this plot
#   geom_point(size = 5) +
#   ggsci::scale_color_npg(name = "Pt. achieved CR") +
#   theme_bw(18) +
#   theme(plot.title = element_text(size = 16, face = "bold")) #-> p.pca.mdanderson.notitle.out

# ggsave("pca_mdanderson_CR_notitle.png", plot = p.pca.mdanderson.notitle.out, width = 7, height = 7, units = "in", dpi = "retina")
# ggsave("pca_mdanderson_CR_notitle.svg", plot = p.pca.mdanderson.notitle.out, width = 7, height = 7, units = "in", dpi = "retina")


#2-2. How about only for mIDh1 patients? 
colData_tp.bl_midh1 <- colData %>% subset(timepoint == "BL") %>% subset(IDH_mutation == "IDH1")
txi_tp.bl_midh1<- tximeta(coldata = colData_tp.bl_midh1)
gse_tp.bl_midh1<- summarizeToGene(txi_tp.bl_midh1)
dds_tp.bl_midh1<- DESeqDataSet(gse_tp.bl_midh1, design = ~achieved_cr)
dds_tp.bl_midh1<- DESeq(dds_tp.bl_midh1)
#vst_tp.bl_midh1<- vst(dds_tp.bl_midh1)
#DESeq2::plotPCA(vst_tp.bl_midh1, intgroup = "achieved_cr", ntop = 100) 

#2-2. How about only for mIDh2 patients? 
colData_tp.bl_midh2 <- colData %>% subset(timepoint == "BL") %>% subset(IDH_mutation == "IDH2")
txi_tp.bl_midh2<- tximeta(coldata = colData_tp.bl_midh2)
gse_tp.bl_midh2<- summarizeToGene(txi_tp.bl_midh2)
dds_tp.bl_midh2<- DESeqDataSet(gse_tp.bl_midh2, design = ~achieved_cr)
dds_tp.bl_midh2<- DESeq(dds_tp.bl_midh2)
#vst_tp.bl_midh2<- vst(dds_tp.bl_midh2)
#DESeq2::plotPCA(vst_tp.bl_midh2, intgroup = "achieved_cr", ntop = 500) 



# RES & gsea --------------------------------------------------------------

#gsea on BL expression profile, annotated with clinical response 
res_cr <- results(dds_tp.bl, contrast = c("achieved_cr", "FALSE", "TRUE")) %>% as.data.frame() 
res_cr$ENSEMBL <- rownames(res_cr) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
map <- AnnotationDbi::select(org.Hs.eg.db, keys = res_cr$ENSEMBL, keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID"))
res_cr <- left_join(res_cr, map, by = "ENSEMBL")

gene.stats_cr <- res_cr$stat
names(gene.stats_cr) <- res_cr$SYMBOL
gene.stats_cr <- gene.stats_cr[!is.na(gene.stats_cr)]

set.seed(69)
tb.gsea_cr <- fgsea::fgseaMultilevel(gs, gene.stats_cr, sampleSize = 500, minSize = 40)
tb.gsea_cr <- tb.gsea_cr[order(tb.gsea_cr$padj), ]

plotEnrichment(gs[["WIERENGA_STAT5A_TARGETS_DN"]], gene.stats_cr) +
  labs(x = "Rank", y = "Enrichment Score",
       title = "WIERENGA_STAT5A_TARGETS_DN",
       subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
  geom_text(aes(10000, -0.15), label = "padj=0.019\nNES=-1.38", size = 10) +
  theme_bw(19) +
  theme(plot.title = element_text(face = "bold", hjust = .5, size = 16), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.WIERENGA_STAT5A_TARGETS_DN.out
ggsave("fgseaMDANDERSON_WIERENGA_STAT5A_TARGETS_DN_CR.png", plot = p.gsea.WIERENGA_STAT5A_TARGETS_DN.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_WIERENGA_STAT5A_TARGETS_DN_CR.svg", plot = p.gsea.WIERENGA_STAT5A_TARGETS_DN.out, height = 5.5, width = 7, units = "in", dpi = "retina")





# plotEnrichment(gs[["WIERENGA_STAT5A_TARGETS_UP"]], gene.stats_cr) +
#   labs(x = "Rank", y = "Enrichment Score",
#        title = "WIERENGA_STAT5A_TARGETS_UP",
#        subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
#   geom_text(aes(29000,.23), label = "padj=0.049\nNES=1.37", size = 10) +
#   theme_bw(19) +
#   theme(plot.title = element_text(face = "bold", hjust = .5, size = 16), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.WIERENGA_STAT5A_TARGETS_UP.out
ggsave("fgseaMDANDERSON_WIERENGA_STAT5A_TARGETS_UP_CR.png", plot = p.gsea.WIERENGA_STAT5A_TARGETS_UP.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_WIERENGA_STAT5A_TARGETS_UP_CR.svg", plot = p.gsea.WIERENGA_STAT5A_TARGETS_UP.out, height = 5.5, width = 7, units = "in", dpi = "retina")

#gsea of patients who achieved CR, comparing REL vs BL
res_tp <- results(dds_cr.true, contrast = c("timepoint", "REL", "BL")) %>% as.data.frame()
res_tp$ENSEMBL <- rownames(res_tp) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist() 
res_tp <- left_join(res_tp, map, by = "ENSEMBL")

gene.stats_tp <- res_tp$stat
names(gene.stats_tp) <- res_tp$SYMBOL
gene.stats_tp <- gene.stats_tp[!is.na(gene.stats_tp)]

set.seed(69)
tb.gsea_tp <- fgsea::fgseaMultilevel(gs, gene.stats_tp, sampleSize = 500, minSize = 40)
tb.gsea_tp <- tb.gsea_tp[order(tb.gsea_tp$padj), ]

plotEnrichment(gs[["STAT5A_04"]], gene.stats_tp) +
  labs(x = "Rank", y = "Enrichment Score",
       title = "Geneset: STAT5A_04",
       subtitle = "Relapse (n=7) vs Baseline (n=11) for patients who achieved CR") +
  geom_text(aes(27000,.23), label = "padj=0.011\nNES=1.37", size = 10) +
  theme_bw(19) +
  theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.STAT5A_04.out
ggsave("fgseaMDANDERSON_STAT5A_04_REL.png", plot = p.gsea.STAT5A_04.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_STAT5A_04_REL.svg", plot = p.gsea.STAT5A_04.out, height = 5.5, width = 7, units = "in", dpi = "retina")


#Responder vs non-responder
res_responder <- results(dds_tp.bl, contrast = c("is_responder", "FALSE", "TRUE")) %>% as.data.frame() 
res_responder$ENSEMBL <- rownames(res_responder) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
map <- AnnotationDbi::select(org.Hs.eg.db, keys = res_responder$ENSEMBL, keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID"))
res_responder <- left_join(res_responder, map, by = "ENSEMBL")

gene.stats_responder <- res_responder$stat
names(gene.stats_responder) <- res_responder$SYMBOL
gene.stats_responder <- gene.stats_responder[!is.na(gene.stats_responder)]

tb.gsea_responder <- fgsea::fgseaMultilevel(gs, gene.stats_responder)
tb.gsea_responder <- tb.gsea_responder[order(tb.gsea_responder$pval), ]




#WRITE FILE
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, "CR_FALSE(n=30)vsTRUE(n=11)")
# openxlsx::addWorksheet(wb, "achivedCR_REL(n=7)_vs_BL(n=11)")
# 
# openxlsx::writeData(wb, "CR_FALSE(n=30)vsTRUE(n=11)", x = res_cr[order(res_cr$padj), ])
# openxlsx::writeData(wb, "achivedCR_REL(n=7)_vs_BL(n=11)", x = res_tp[order(res_tp$padj), ])
# 
# openxlsx::saveWorkbook(wb, "MDAnderson_PtData_DESeq2_results.csv")
# GSEA HOX GENES ----------------------------------------------------------

#CR 
set.seed(69)
tb.gsea.hox_cr <- fgsea::fgseaMultilevel(gs.hox, gene.stats_cr, sampleSize = 500, minSize = 40)
tb.gsea.hox_cr <- tb.gsea.hox_cr[order(tb.gsea.hox_cr$padj), ]


# plotEnrichment(gs.hox[["CHEN_HOXA5_TARGETS_9HR_UP"]], gene.stats_cr) + 
#   labs(x = "Rank", y = "Enrichment Score", 
#        title = "Geneset: CHEN_HOXA5_TARGETS_9HR_UP", 
#        subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
#   geom_text(aes(30000, 0.4), label = "padj<0.01\nNES=2.24", size = 10) +
#   theme_bw(19) +
#   theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_CR.out
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_UP_CR.png", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_CR.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_UP_CR.svg", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_CR.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")

# plotEnrichment(gs.hox[["CHEN_HOXA5_TARGETS_9HR_DN"]], gene.stats_cr) + 
#   labs(x = "Rank", y = "Enrichment Score", 
#        title = "Geneset: CHEN_HOXA5_TARGETS_9HR_UP", 
#        subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
#   geom_text(aes(10000, -.4), label = "padj<0.01\nNES=-1.92", size = 10) +
#   theme_bw(19) +
#   theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_CR.out
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_DN_CR.png", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_CR.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_DN_CR.svg", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_CR.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")

plotEnrichment(gs.hox[["TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP"]], gene.stats_cr) +
  labs(x = "Rank", y = "Enrichment Score",
       title = "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP",
       subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
  geom_text(aes(30000, .4), label = "padj=1.85E-9\nNES=2.15", size = 10) +
  theme_bw(19) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = .5), plot.subtitle = element_text(size = 16, hjust = .5)) -> p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.out
ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.png", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.svg", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.out, height = 5.5, width = 7, units = "in", dpi = "retina")


plotEnrichment(gs.hox[["TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP"]], gene.stats_cr) +
  labs(x = "Rank", y = "Enrichment Score",
       title = "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP",
       subtitle = "Pt. did not achieve CR (n=30) vs Pt. achieved CR (n=11)") +
  geom_text(aes(30000, .4), label = "padj=3.18E-6\nNES=2.31", size = 10) +
  theme_bw(19) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = .5), plot.subtitle = element_text(size = 16, hjust = .5)) -> p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP_CR.out
ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP_CR.png", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP_CR.svg", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_CR.out, height = 5.5, width = 7, units = "in", dpi = "retina")

##REL
set.seed(69)
tb.gsea.hox_tp <- fgsea::fgseaMultilevel(gs.hox, gene.stats_tp, sampleSize = 500, minSize = 40)
tb.gsea.hox_tp <- tb.gsea.hox_tp[order(tb.gsea.hox_tp$padj), ]

# plotEnrichment(gs.hox[["CHEN_HOXA5_TARGETS_9HR_UP"]], gene.stats_tp) + 
#   labs(x = "Rank", y = "Enrichment Score", 
#        title = "Geneset: CHEN_HOXA5_TARGETS_9HR_UP", 
#        subtitle = "Pt. achieved CR: RELAPSE(n=7) vs BASELINE(n=11)") +
#   geom_text(aes(27000, 0.4), label = "padj<0.001\nNES=2.52", size = 10) +
#   theme_bw(19) +
#   theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_REL.out
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_UP_REL.png", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_REL.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_UP_REL.svg", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_UP_REL.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")

# plotEnrichment(gs.hox[["CHEN_HOXA5_TARGETS_9HR_DN"]], gene.stats_tp) + 
#   labs(x = "Rank", y = "Enrichment Score", 
#        title = "Geneset: CHEN_HOXA5_TARGETS_9HR_DN", 
#        subtitle = "Pt. achieved CR: RELAPSE(n=7) vs BASELINE(n=11)") +
#   geom_text(aes(10000, -0.3), label = "padj=0.055\nNES=-1.48", size = 10) +
#   theme_bw(19) +
#   theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 16)) -> p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_REL.out
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_DN_REL.png", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_REL.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")
# ggsave("fgseaMDANDERSON_CHEN_HOXA5_TARGETS_9HR_DN_REL.svg", plot = p.gsea.CHEN_HOXA5_TARGETS_9HR_DN_REL.out, height = 5.5, width = 7.5, units = "in", dpi = "retina")

plotEnrichment(gs.hox[["TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP"]], gene.stats_tp) +
  labs(x = "Rank", y = "Enrichment Score",
       title = "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP",
       subtitle = "Relapse (n=7) vs Baseline (n=11) for patients who achieved CR") +
  geom_text(aes(27000, .3), label = "padj=3.04E-6\nNES=1.87", size = 10) +
  theme_bw(19) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = .5), plot.subtitle = element_text(size = 16, hjust = .5)) -> p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_REL.out

ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_REL.png", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_REL.out, height = 5.5, width = 7, units = "in", dpi = "retina")
ggsave("fgseaMDANDERSON_TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_REL.svg", plot = p.gsea.TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_UP_REL.out, height = 5.5, width = 7, units = "in", dpi = "retina")


##RESPONDERS
tb.gsea.hox_responder <- fgsea::fgseaMultilevel(gs.hox, gene.stats_responder, sampleSize = 500, minSize = 40)
tb.gsea.hox_responder <- tb.gsea.hox_responder[order(tb.gsea.hox_responder$padj), ]


###Paired samples
res_tp_paired <- results(dds.paired, contrast = c("timepoint", "REL", "BL")) %>% as.data.frame
res_tp_paired$ENSEMBL <- rownames(res_tp_paired) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist() 
res_tp_paired <- left_join(res_tp_paired, map, by = "ENSEMBL")



###HOX genes exploration


map[grep("HOX", map$SYMBOL),]$SYMBOL

x <- DESeq2::counts(dds_tp.bl, normalized = TRUE) %>% as.data.frame()
x$ENSEMBL <- rownames(x) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
x <- left_join(x, map)
x <- x[x$SYMBOL %in% map[grep("HOX", map$SYMBOL),]$SYMBOL, ]  #Select hox genes
x.meta <- colData_tp.bl %>% dplyr::select(names, achieved_cr)
x.meta <- x.meta[order(x.meta$achieved_cr),]

x %>% 
  pivot_longer(cols = 1:41, names_to = "names", values_to = "deseq_count") %>% 
  left_join(x.meta) -> x


int_symbols <- x %>% group_by(SYMBOL) %>% dplyr::summarise(sum = sum(deseq_count)) %>%
  .[order(.$sum), ] %>% 
  subset(sum >10) %>%
  .$SYMBOL
x %>% 
  subset(SYMBOL %in% int_symbols) %>% 
  ggplot(aes(achieved_cr, deseq_count, fill = achieved_cr)) +
  geom_boxplot() +
  facet_wrap(~SYMBOL, scales = "free_y")


# Plot counts -------------------------------------------------------------

counts_tp.bl <- counts(dds_tp.bl, normalized = TRUE)
counts_cr.true <- counts(dds_cr.true, normalized = TRUE)

deseq_counts_ensembl_to_symbol <- function(deseq_count){
  deseq_count <- as.data.frame(deseq_count)
  deseq_count$ENSEMBL <- rownames(deseq_count) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
  deseq_count <- left_join(deseq_count, map)
  return(deseq_count)
}

counts_tp.bl <- deseq_counts_ensembl_to_symbol(counts_tp.bl)
counts_cr.true <- deseq_counts_ensembl_to_symbol(counts_cr.true)


# -------------------------------------------------------------------------

#LSC17
# counts[SYMBOL == "DNMT3B", 1:51] * 0.0874 +
#   counts[SYMBOL == "ZBTB46", 1:51] * 0.0874 +
gene_list <- c("DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B", 
               "MMRN1", "DPYSL3", "FAM30A", "CDK6", "CPXM1", 
               "SOCS2", "SMIM24", "EMP1", "BEX3", "CD34", "AKR1C3", "ADGRG1")

score_list <- c(.0874, -.0347, .00865, -.0138, .00582, .0258, .0284, .0196,
                -.0704, -.0258, .0271, -.0226, .0146, .0465, .0338, -.0402, .0501)

score <- vector(mode = "numeric", length = 51)
for(i in 1:length(gene_list)){
  score <- score + (as.numeric(counts[SYMBOL == gene_list[i], 1:51]) * score_list[i])
}
  
  
dt <- data.table(colData)
dt$lsc17 <- score

dt$absolute_cr <- factor((dt$Best_Response == "CR"), levels = c(TRUE, FALSE))

glm(achieved_cr ~ lsc17, data = dt, family = "binomial") %>% summary()
glm(absolute_cr ~ lsc17, data = dt, family = "binomial") %>% summary()


dt$hoxa2 <- counts[SYMBOL == "", 1:51] %>% as.numeric() 
glm(achieved_cr ~ hoxa2, data = dt, family = "binomial") %>% summary()
glm(absolute_cr ~ hoxa2, data = dt, family = "binomial") %>% summary()




# MDM2 --------------------------------------------------------------------


rownames(counts) <- rownames(counts) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
df.mdm2 <- counts %>% subset(rownames(.) == "ENSG00000135679") %>% t() %>% cbind(colData, .)


results(dds.bl, contrast = c("is_responder", "FALSE", "TRUE")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000135679")


results(dds, contrast = c("timepoint", "REL", "BL")) %>% 
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000135679")

df.mdm2 %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(Best_Response, ENSG00000135679)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Best Response", y = "MDM2 - DESeq2 Normalized Counts",
       title = "Counts vs Best Response") +
  theme_bw() -> p1

df.mdm2 %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(is_responder, ENSG00000135679)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Responder", y = NULL,
       title = "Nonresponder vs responder") +
  geom_label(data = NULL, aes(x = 1.5, y = 11000), label = "L2FC = .304\npval = .083\npadj=.72") +
  theme_bw()-> p2

df.mdm2 %>% 
  ggplot(aes(timepoint, ENSG00000135679)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Timepoint", y = NULL,
       title = "REL vs BL") +
  geom_label(data = NULL, aes(x = 1.5, y = 11000), label = "L2FC = .198\npval = .344\npadj=.86") +
  theme_bw() -> p3



res.paired.rel %>% subset(rownames(.) == "ENSG00000135679")

df.mdm2.paired <- df.mdm2 %>% subset(id %in% id.paired) 
  
ggplot(df.mdm2.paired, aes(timepoint, ENSG00000135679, shape = id)) +
  geom_point() +
  scale_shape_manual(values = c(1:7)) +
  labs(x = "Timepoint", y = NULL,
       title = "Paired samples only") +
  geom_label(data = NULL, aes(x = 1.5, y = 8000), label = "L2FC = .3\npval = .33\npadj=.83") +
  theme_bw() +
  geom_path(aes(group = id)) -> p4

(p1 + p2) + (p3 + p4) -> p.out

ggsave(filename = "MDM2_email.png",
       plot = p.out, width = 18, height = 6, dpi = "retina", units = "in")



# check socs2 --------------------------------------------------------------

rownames(counts) <- rownames(counts) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()
df.socs2 <- counts %>% subset(rownames(.) == "ENSG00000120833") %>% t() %>% cbind(colData, .)


results(dds.bl, contrast = c("is_responder", "FALSE", "TRUE")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000120833")


results(dds, contrast = c("timepoint", "REL", "BL")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000120833")

results(dds.paired, contrast = c("timepoint", "REL", "BL")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000120833")

df.socs2 %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(Best_Response, ENSG00000120833)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Best Response", y = "SOCS2 - DESeq2 Normalized Counts",
       title = "Counts vs Best Response") +
  theme_bw() -> p.socs2.1

df.socs2 %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(is_responder, ENSG00000120833)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  geom_label(data = NULL, aes(1.5, 2000, label = "l2fc = .523\npval=.145\npadj=.787")) +
  labs(x = "Responder", y = NULL,
       title = "Nonresponder vs responder") +
  theme_bw() -> p.socs2.2

df.socs2 %>% 
  ggplot(aes(timepoint, ENSG00000120833)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  geom_label(data = NULL, aes(1.5, 2000, label = "l2fc = -.251\npval=.614\npadj=.941")) +
  labs(x = "Timepoint", y = NULL,
       title = "REL vs BL") +
  theme_bw() -> p.socs2.3


#
df.socs2.paired <- df.socs2 %>% subset(id %in% id.paired) 

ggplot(df.socs2.paired, aes(timepoint, ENSG00000120833, shape = id)) +
  geom_point() +
  scale_shape_manual(values = c(1:7)) +
  labs(x = "Timepoint", y = NULL,
       title = "Paired samples only") +
  geom_label(data = NULL, aes(1.5, 600, label = "l2fc = .32\npval=.607\npadj=.925")) +
  theme_bw() +
  geom_path(aes(group = id)) -> p.socs2.4


(p.socs2.1 + p.socs2.2) + (p.socs2.3 + p.socs2.4) -> p.socs2.out
ggsave(filename = "socs2_email.png",
       plot = p.socs2.out, width = 18, height = 6, dpi = "retina", units = "in")

# check il2ra -------------------------------------------------------------
df.il2ra <- counts %>% subset(rownames(.) == "ENSG00000134460") %>% t() %>% cbind(colData, .)


#
results(dds.bl, contrast = c("is_responder", "FALSE", "TRUE")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000134460")


results(dds, contrast = c("timepoint", "REL", "BL")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000134460")

results(dds.paired, contrast = c("timepoint", "REL", "BL")) %>%
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.)) %>% 
  mutate(ENSEMBL = ENSEMBL %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  subset(ENSEMBL == "ENSG00000134460")
#

df.il2ra %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(Best_Response, ENSG00000134460)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Best Response", y = "IL2RA - DESeq2 Normalized Counts",
       title = "Counts vs Best Response") +
  theme_bw() -> p.il2ra.1


df.il2ra %>%
  subset(timepoint == "BL") %>%
  ggplot(aes(is_responder, ENSG00000134460)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Responder", y = NULL,
       title = "Nonresponder vs responder") +
  geom_label(data = NULL, aes(1.5, 3000, label = "l2fc = -.523\npval=.363\npadj=.882")) +
  theme_bw() -> p.il2ra.2

df.il2ra %>% 
  ggplot(aes(timepoint, ENSG00000134460)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = .05) +
  labs(x = "Timepoint", y = NULL,
       title = "REL vs BL") +
  geom_label(data = NULL, aes(1.5, 3000, label = "l2fc = -2.36\npval=<0.001\npadj=.187")) +
  theme_bw() #-> p.il2ra.3

#
df.il2ra.paired <- df.il2ra %>% subset(id %in% id.paired) 

ggplot(df.il2ra.paired, aes(timepoint, ENSG00000134460, shape = id)) +
  geom_point() +
  scale_shape_manual(values = c(1:7)) +
  labs(x = "Timepoint", y = NULL,
       title = "Paired samples only") +
  theme_bw() +
  geom_label(data = NULL, aes(2, 3300, label = "l2fc = -2.78\npval=<0.001\npadj=.022")) +
  geom_path(aes(group = id)) -> p.il2ra.4

(p.il2ra.1 + p.il2ra.2) + (p.il2ra.3 + p.il2ra.4) -> p.il2ra.out

# ggsave(filename = "il2ra_email.png",
#        plot = p.il2ra.out, width = 18, height = 6, dpi = "retina", units = "in")



# manuscript il2ra plot ---------------------------------------------------
library(ggtext)
df.il2ra %>% 
  ggplot(aes(timepoint, ENSG00000134460, fill = timepoint)) +
  #geom_boxplot(outlier.shape = NA, alpha = .4) + 
  geom_violin(alpha = .4) +
  geom_jitter(width = .05) +
  labs(x = "Timepoint", y = NULL,
       title = "IL2RA - DESeq2 normalized counts") +
  geom_richtext(data = NULL, 
                aes(1.5, 3000, 
                    label = "l2fc = -2.36<br>
                    pval = <0.001<br>padj = 0.187"),
                size = 6, hjust = 0,
                color = "navy", fill = "white"
                ) +
  scale_fill_manual(values = c(BL = "navy", REL = "firebrick"),
                    guide = "none") +
  theme_classic(13) -> p1


#
df.il2ra.paired <- df.il2ra %>% subset(id %in% id.paired) 

ggplot(df.il2ra.paired, aes(timepoint, ENSG00000134460, shape = id)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1:7), name = "Patient id") +
  labs(x = "Timepoint", y = NULL,
       title = "Paired patient samples only") +
  theme_classic(13) +
  geom_richtext(data = NULL, 
                aes(1.5, 3300, label = "l2fc = -2.78<br>
                    pval = <0.001<br>padj = 0.022"),
                color = "navy", fill = "white", size = 6, hjust = 0) +
  geom_path(aes(group = id)) -> p2

p1 + p2 + plot_layout(widths = c(1, 1.2)) -> p.out

ggsave(filename = "Il2ra_count.png", plot = p.out, 
       height = 4.43, width = 10, units = "in", dpi = 300)
ggsave(filename = "Il2ra_count.svg", plot = p.out, 
       height = 4.43, width = 10, units = "in", dpi = 300)



# -------------------------------------------------------------------------
#2021.11.04 Exploring count data
#
library(data.table)
natcom_counts <- counts(dds, normalized = TRUE)
natcom_counts <- natcom_counts %>% 
  as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) 

natcom_counts %>% 
  data.table() %>% 
  pivot_longer(starts_with("SRR"), names_to = "names", values_to = "deseq2_counts") %>% 
  left_join(colData) -> natcom_counts

AnnotationDbi::select(org.Hs.eg.db,
                      keys = c("CDKN3"), 
                      keytype = "SYMBOL",
                      columns = "ENSEMBL") -> temp

natcom_counts %>% 
  subset(ENSEMBL %in% "ENSG00000260314") %>% 
  ggplot(aes(achieved_cr, deseq2_counts)) + 
  geom_boxplot() +
  facet_wrap(~ENSEMBL, scales = "free_y")


results(dds, contrast = c("achieved_cr", "FALSE", "TRUE")) %>% as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  left_join(map) %>% 
  .[order(.$pval), ] %>% 
  subset(baseMean > 1000)

results(dds.paired, contrast = c("timepoint", "REL", "BL")) %>% as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  left_join(map) %>% 
  .[order(.$pval), ] %>% 
  subset(baseMean > 1000)

results(dds.bl, contrast = c("is_responder", "FALSE", "TRUE")) %>% as.data.frame() %>% 
  mutate(ENSEMBL = rownames(.) %>% str_split(pattern = "\\.") %>% lapply(function(x){x[1]}) %>% unlist()) %>% 
  left_join(map) %>% 
  .[order(.$pval), ] %>% 
  subset(baseMean > 1000)
