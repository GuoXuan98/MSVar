# Main analysis on HCC tumor samples.

library(MSVar)
library(openxlsx)
library(DMwR2)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(pheatmap)
library(clusterProfiler)

# Read data ---------------------------------------------------------------

# Raw intensity and batch information for tumor samples.
HCC.raw.intensity.tumor <- read.xlsx("data/HCC_159Tumor.xlsx",
                                     rowNames = T, check.names = F, sheet = 1)
HCC.batch.info.tumor <- read.xlsx("data/HCC_159Tumor.xlsx",
                                  rowNames = F, check.names = F, sheet = 2)

# Raw intensity and batch information for normal samples.
HCC.raw.intensity.normal <- read.xlsx("data/HCC_159Normal.xlsx",
                                      rowNames = T, check.names = F, sheet = 1)
HCC.batch.info.normal <- read.xlsx("data/HCC_159Normal.xlsx",
                                   rowNames = F, check.names = F, sheet = 2)

# Clinical information for 159 HCC patients.
HCC.clinical.info <- read.xlsx("data/HCC_159Tumor.xlsx",
                               rowNames = T, check.names = F, sheet = 3)

# HCC, tumor -------------------------------------------------------------

# Construct a proObj for the group of tumor samples.
HCC.tumor <- proObjFromTMT(HCC.raw.intensity.tumor, HCC.batch.info.tumor,
                           IDs = rownames(HCC.raw.intensity.tumor),
                           name = "HCC.tumor")

# Perform technical variance estimation procedure.
HCC.tumor <- estTechVar(HCC.tumor)
# Perform biological variance estimation procedure.
HCC.tumor <- estBioVar(HCC.tumor)

# Derive posterior M-value.
HCC.tumor <- PostM(HCC.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across HCC tumor samples.
HCC.HVP <- VarTest(HCC.tumor)


# HCC, normal -------------------------------------------------------------

# Construct a proObj for the group of normal samples.
HCC.normal <- proObjFromTMT(HCC.raw.intensity.normal, HCC.batch.info.normal,
                            IDs = rownames(HCC.raw.intensity.normal),
                            name = "HCC.normal")

# Perform technical variance estimation procedure.
HCC.normal <- estTechVar(HCC.normal)
# Perform biological variance estimation procedure.
HCC.normal <- estBioVar(HCC.normal)

# Derive posterior M-value.
HCC.normal <- PostM(HCC.normal)


# DVP identification ------------------------------------------------------

# Conduct hypothesis testing procedures for identifying tumor-high DVPs.
HCC.tumor.DVP.test <- DiffVar(HCC.tumor, HCC.normal)

# Identify DVPs under the condition that the BH-corrected P-value is < 0.01.
HCC.tumor.DVP <- rownames(HCC.tumor.DVP.test)[which(HCC.tumor.DVP.test[, "padj"] < 0.01)]


# Expression extraction ---------------------------------------------------

# Posterior M-value.
PostM.HCC.tumor <- HCC.tumor$PostM
# Z-score scaled posterior M-value.
PostZScr.HCC.tumor <- HCC.tumor$PostZScr

# Imputation --------------------------------------------------------------

# Proteins whose biological variance estimates are clearly larger than 0.
UpperPart <- attr(HCC.tumor$resBio$est.bio.res, "group1")

# Imputation among proteins in the upper part.
set.seed(110)
PostZScrUpperImp.HCC.tumor <- knnImputation(PostZScr.HCC.tumor[UpperPart, ],
                                            k = 10, scale = F)

# Imputation among all proteins.
set.seed(110)
PostZScrImp.HCC.tumor <- knnImputation(PostZScr.HCC.tumor,
                                       k = 10, scale = F)

# Sample classification ---------------------------------------------------

# Consensus clustering by 1747 significant tumor-high DVPs.
HCC.tumor.DVP.CCRes <- ConsensusClusterPlus(
  as.matrix(PostZScrUpperImp.HCC.tumor[HCC.tumor.DVP, ]),
  maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "km", title = "./CSSC_res/",
  distance = "euclidean", seed = 110, plot = "png")

# HCC tumor samples are clustered into 2 groups: Group1 & Group2.
HCC.tumor.DVP.2CRes <- HCC.tumor.DVP.CCRes[[2]]$consensusClass

# Extract samples in each group.
HCC.G1 <- names(HCC.tumor.DVP.2CRes)[HCC.tumor.DVP.2CRes == 1]
HCC.G2 <- names(HCC.tumor.DVP.2CRes)[HCC.tumor.DVP.2CRes == 2]


# DEP identification ------------------------------------------------------

# MSVar fits the single-group model separately for Group1 & Group2, and then
# conduct differential expression analysis between Group1 & Group2.

# Extract batch information and intensity for samples in Group1.
HCC.batch.info.G1 <- t(sapply(1:nrow(HCC.batch.info.tumor), function(i) {
  temp <- HCC.batch.info.tumor[i, ] %in% HCC.G1
  temp[1] <- T

  sample.info <- HCC.batch.info.tumor[i, ]
  sample.info[!temp] <- ""

  unlist(sample.info)
}))
HCC.raw.intensity.G1 <- HCC.raw.intensity.tumor[, as.vector(t(HCC.batch.info.G1))[nchar(t(HCC.batch.info.G1)) > 0]]

# Extract batch information and intensity for samples in Group2.
HCC.batch.info.G2 <- t(sapply(1:nrow(HCC.batch.info.tumor), function(i) {
  temp <- HCC.batch.info.tumor[i, ] %in% HCC.G2
  temp[1] <- T

  sample.info <- HCC.batch.info.tumor[i, ]
  sample.info[!temp] <- ""

  unlist(sample.info)
}))
HCC.raw.intensity.G2 <- HCC.raw.intensity.tumor[, as.vector(t(HCC.batch.info.G2))[nchar(t(HCC.batch.info.G2)) > 0]]

# Construct a proObj for the Group1 of tumor samples.
HCC.tumor.G1 <- proObjFromTMT(HCC.raw.intensity.G1, HCC.batch.info.G1,
                              IDs = rownames(HCC.raw.intensity.G1),
                              name = "HCC.tumor.G1")

# Perform technical & biological variance estimation procedure.
HCC.tumor.G1 <- estTechVar(HCC.tumor.G1)
HCC.tumor.G1 <- estBioVar(HCC.tumor.G1)

# Construct a proObj for the Group1 of tumor samples.
HCC.tumor.G2 <- proObjFromTMT(HCC.raw.intensity.G2, HCC.batch.info.G2,
                              IDs = rownames(HCC.raw.intensity.G2),
                              name = "HCC.tumor.G2")

# Perform technical & biological variance estimation procedure.
HCC.tumor.G2 <- estTechVar(HCC.tumor.G2)
HCC.tumor.G2 <- estBioVar(HCC.tumor.G2)

# Conduct hypothesis testing procedures for identifying DEPs between Group1 & Group2.
DEP.G1High.Res <- DiffMean(HCC.tumor.G1, HCC.tumor.G2, scaleV.raw.FC = 0,
                           alternative = "upper", p.adjust.method = "BH")
DEP.G2High.Res <- DiffMean(HCC.tumor.G2, HCC.tumor.G1, scaleV.raw.FC = 0,
                           alternative = "upper", p.adjust.method = "BH")

# Identify DEPs under the condition that the BH-corrected P-value is < 0.01.
DEP.G1High <- rownames(DEP.G1High.Res)[which(DEP.G1High.Res$padj < 0.01)]
DEP.G1High <- intersect(rownames(DEP.G1High.Res[order(DEP.G1High.Res$padj), ]),
                        DEP.G1High)
DEP.G2High <- rownames(DEP.G2High.Res)[which(DEP.G2High.Res$padj < 0.01)]
DEP.G2High <- intersect(rownames(DEP.G2High.Res[order(DEP.G2High.Res$padj, decreasing = T), ]),
                        DEP.G2High)


# Survival analysis -------------------------------------------------------

# Extract overall survival information of HCC patients.
HCC.surv.OS <- data.frame(time = HCC.clinical.info[, "OS_days"],
                          status = HCC.clinical.info[, "OS_event"],
                          row.names = rownames(HCC.clinical.info))
HCC.surv.OS$subgroup <- HCC.tumor.DVP.2CRes[rownames(HCC.surv.OS)]

# Create survival curve.
HCC.surv.fit.Res <- survfit(Surv(time, status) ~ subgroup,
                            data = HCC.surv.OS)

# Compute p-value by using the log-rank test.
HCC.surv.diff.Pval <- signif(surv_pvalue(HCC.surv.fit.Res, data = HCC.surv.OS)$pval, 3)
Leg <- sapply(1:2, function(i)
  paste0("G", i, " (", length(which(HCC.surv.OS$subgroup == i)),
         " patients)"))

# Plot Kaplan-Meier survival curve.
HCC.sur.curve <- ggsurvplot(
  fit = HCC.surv.fit.Res, data = HCC.surv.OS,
  pval = HCC.surv.diff.Pval,
  pval.coord = c(50*0.6, 0.1),
  pval.size = 8, size = 2.2,
  xlab = "Overall survival (month)", ylab = "Survival ratio",
  title = "1747 DVPs (BH-adjusted p < 0.01)",
  xlim = c(0, 50),
  legend = c(0.3, 0.33), legend.title = "", legend.labs = Leg,
  palette = c("#FED439FF", "#709AE1FF"),
  font.x = 25, font.y = 25,
  font.legend = 20, font.tickslab = 20,
  ggtheme = theme(aspect.ratio = 1/1,
                  panel.grid = element_blank(), panel.background = element_blank(),
                  axis.ticks = element_blank(), axis.text = element_blank(),
                  axis.ticks.length.x = unit(0.3, 'cm'),
                  axis.ticks.length.y = unit(0.3, 'cm'),
                  axis.title.x = element_text(hjust = 0.5, vjust = -0.3, colour = "black"),
                  axis.text.x = element_text(hjust = 0.5, vjust = -0.2, colour = "black"),
                  axis.title.y = element_text(angle = 90,
                                              hjust = 0.5, vjust = 2, colour = "black"),
                  axis.text.y = element_text(angle = 90,
                                             hjust = 0.5, vjust = 1.2, colour = "black"),
                  plot.margin = margin(t = 14, r = 9, b = 15, l = 8),
                  plot.title = element_text(size = 26),
                  legend.spacing.y = unit(1, 'cm'),
                  panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2)))


# Heatmap -----------------------------------------------------------------

# Obtain the consensus clustering tree.
clus.order <- HCC.tumor.DVP.CCRes[[2]]$consensusTree$order
bk <- seq(-0.6, 0.6, by = 0.05)

HCC.clinical.category <- HCC.clinical.info[, 1:4]
HCC.clinical.category$`BCLC.stage` <- as.character(HCC.clinical.info$BCLC.stage)
HCC.clinical.category$`AFP.more200` <- HCC.clinical.info$`Preoperative.AFP（ng/mL）` > 200

# Specify the annotations of columns.
anno.col.clinic <- data.frame("AFP more200" = as.factor(HCC.clinical.category$`AFP.more200`),
                              "BCLC stage" = as.character(HCC.clinical.category$`BCLC.stage`),
                              ConClus = as.factor(HCC.tumor.DVP.2CRes),
                              row.names = names(HCC.tumor.DVP.2CRes),
                              check.names = F)
anno.col.clinic <- anno.col.clinic[clus.order, ]
anno.col.clinic <- anno.col.clinic[order(anno.col.clinic$ConClus), ]

# Specify annotation_col track colors.
anno.color.clinic <- list(
  "AFP more200" = c("TRUE" = "black", "FALSE" = "gray"),
  "BCLC stage" = c("A" = "pink", "B" = "tomato", "C" = "brown"),
  ConClus = c("1" = "#FED439FF", "2" = "#709AE1FF"))

# Plot heatmap of DEPs between 2 subgroups.
HCC.DEP.heatmap <- pheatmap(
  PostZScrImp.HCC.tumor[c(DEP.G1High, DEP.G2High), rownames(anno.col.clinic)],
  show_rownames = F, show_colnames = F,
  cluster_cols = F, cluster_row = F,
  gaps_col = length(HCC.G1), gaps_row = length(DEP.G1High),
  main = paste0("G1 (", length(HCC.G1), " patients) & ",
                "G2 (", length(HCC.G2), " patients)"),
  annotation_col = anno.col.clinic, annotation_colors = anno.color.clinic,
  col = c(colorRampPalette(colors = c("blue", "white"))(round(length(bk)/2)),
          colorRampPalette(colors = c("white", "red"))(ceiling(length(bk)/2))),
  breaks = bk)


# Enrichment analysis -----------------------------------------------------

# KEGG enrichment analysis of Group1-high DEPs.
DEP.G1High.ID <- bitr(DEP.G1High,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

DEP.G1High.Enrich.Kegg <- enrichKEGG(gene = DEP.G1High.ID$ENTREZID,
                                     pvalueCutoff = 0.5, qvalueCutoff = 0.5)

# KEGG enrichment analysis of Group2-high DEPs.
DEP.G2High.ID <- bitr(DEP.G2High,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

DEP.G2High.Enrich.Kegg <- enrichKEGG(gene = DEP.G2High.ID$ENTREZID,
                                     pvalueCutoff = 0.5, qvalueCutoff = 0.5)

# Feature association evaluating ------------------------------------------

# Evaluating the association between protein expression and patient prognosis.

# Fit a Cox proportional hazards regression model for each protein,
# and the two-tailed p-value produced by it serves as the statistical measure of association significance.
PostM.HCC.tumor.Cox.Res <- exprCox(
  expr = PostM.HCC.tumor[, rownames(HCC.clinical.info)],
  time = HCC.clinical.info$OS_days,
  status = HCC.clinical.info$OS_event)

# Evaluating the association between protein expression and BCLC stage.

# Conduct Welch’s one-way test or ANOVA for categorical feature,
# and the resulting p-value serves as the statistical measure of protein-feature association.
PostM.HCC.tumor.BCLC.Res <- exprCatFeature(expr = PostM.HCC.tumor[, rownames(HCC.clinical.category)],
                                           feature = HCC.clinical.category[, "BCLC.stage"])

# Evaluating the association between protein expression and indicator of abnormally high AFP level.
PostM.HCC.tumor.AFP.Res <- exprCatFeature(expr = PostM.HCC.tumor[, rownames(HCC.clinical.category)],
                                          feature = HCC.clinical.category[, "AFP.more200"])
