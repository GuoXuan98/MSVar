# Main analysis on LUAD tumor samples.

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
LUAD.raw.intensity.tumor <- read.xlsx("data/LUAD_110Tumor.xlsx",
                                      rowNames = T, check.names = F, sheet = 1)
LUAD.batch.info.tumor <- read.xlsx("data/LUAD_110Tumor.xlsx",
                                   rowNames = F, check.names = F, sheet = 2)

# Raw intensity and batch information for normal samples.
LUAD.raw.intensity.normal <- read.xlsx("data/LUAD_101Normal.xlsx",
                                       rowNames = T, check.names = F, sheet = 1)
LUAD.batch.info.normal <- read.xlsx("data/LUAD_101Normal.xlsx",
                                    rowNames = F, check.names = F, sheet = 2)

# Clinical information for 110 LUAD patients.
LUAD.clinical.info <- read.xlsx("data/LUAD_110Tumor.xlsx",
                                rowNames = T, check.names = F, sheet = 3)

# LUAD, tumor -------------------------------------------------------------

# Construct a proObj for the group of tumor samples.
LUAD.tumor <- proObjFromTMT(LUAD.raw.intensity.tumor, LUAD.batch.info.tumor,
                            IDs = rownames(LUAD.raw.intensity.tumor),
                            name = "LUAD.tumor")

# Perform technical variance estimation procedure.
LUAD.tumor <- estTechVar(LUAD.tumor)
# Perform biological variance estimation procedure.
LUAD.tumor <- estBioVar(LUAD.tumor)

# Derive posterior M-value.
LUAD.tumor <- PostM(LUAD.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across LUAD tumor samples.
LUAD.HVP <- VarTest(LUAD.tumor)


# LUAD, normal -------------------------------------------------------------

# Construct a proObj for the group of normal samples.
LUAD.normal <- proObjFromTMT(LUAD.raw.intensity.normal, LUAD.batch.info.normal,
                             IDs = rownames(LUAD.raw.intensity.normal),
                             name = "LUAD.normal")

# Perform technical variance estimation procedure.
LUAD.normal <- estTechVar(LUAD.normal)
# Perform biological variance estimation procedure.
LUAD.normal <- estBioVar(LUAD.normal)

# Derive posterior M-value.
LUAD.normal <- PostM(LUAD.normal)


# DVP identification ------------------------------------------------------

# Conduct hypothesis testing procedures for identifying tumor-high DVPs.
LUAD.tumor.DVP.test <- DiffVar(LUAD.tumor, LUAD.normal)

# Identify DVPs under the condition that the BH-corrected P-value is < 0.01.
LUAD.tumor.DVP <- rownames(LUAD.tumor.DVP.test)[which(LUAD.tumor.DVP.test[, "padj"] < 0.01)]


# Expression extraction ---------------------------------------------------

# Posterior M-value.
PostM.LUAD.tumor <- LUAD.tumor$PostM
# Z-score scaled posterior M-value.
PostZScr.LUAD.tumor <- LUAD.tumor$PostZScr

# Imputation --------------------------------------------------------------

# Proteins whose biological variance estimates are clearly larger than 0.
UpperPart <- attr(LUAD.tumor$resBio$est.bio.res, "group1")

# Imputation among proteins in the upper part.
set.seed(110)
PostZScrUpperImp.LUAD.tumor <- knnImputation(PostZScr.LUAD.tumor[UpperPart, ],
                                             k = 10, scale = F)
# Imputation among all proteins.
set.seed(110)
PostZScrImp.LUAD.tumor <- knnImputation(PostZScr.LUAD.tumor,
                                        k = 10, scale = F)

# Sample classification ---------------------------------------------------

# Consensus clustering by 740 significant tumor-high DVPs.
LUAD.tumor.DVP.CCRes <- ConsensusClusterPlus(
  as.matrix(PostZScrUpperImp.LUAD.tumor[LUAD.tumor.DVP, ]),
  maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "km", title = "./CSSC_res/",
  distance = "euclidean", seed = 110, plot = "png")

# LUAD tumor samples are clustered into 3 groups: Group1, Group2 & Group3.
LUAD.tumor.DVP.3CRes <- LUAD.tumor.DVP.CCRes[[3]]$consensusClass

# Extract samples in each group.
LUAD.G1 <- names(LUAD.tumor.DVP.3CRes)[LUAD.tumor.DVP.3CRes == 1]
LUAD.G2 <- names(LUAD.tumor.DVP.3CRes)[LUAD.tumor.DVP.3CRes == 2]
LUAD.G3 <- names(LUAD.tumor.DVP.3CRes)[LUAD.tumor.DVP.3CRes == 3]


# DEP identification ------------------------------------------------------

# MSVar fits the single-group model separately for Group1 & Group2, and then
# conduct differential expression analysis between Group1 & Group2.

# Extract batch information and intensity for samples in Group1.
LUAD.batch.info.G1 <- t(sapply(1:nrow(LUAD.batch.info.tumor), function(i) {
  temp <- LUAD.batch.info.tumor[i, ] %in% LUAD.G1
  temp[6] <- T

  sample.info <- LUAD.batch.info.tumor[i, ]
  sample.info[!temp] <- ""

  unlist(sample.info)
}))
LUAD.raw.intensity.G1 <- LUAD.raw.intensity.tumor[, as.vector(t(LUAD.batch.info.G1))[nchar(t(LUAD.batch.info.G1)) > 0]]

# Extract batch information and intensity for samples in Group2.
LUAD.batch.info.G2 <- t(sapply(1:nrow(LUAD.batch.info.tumor), function(i) {
  temp <- LUAD.batch.info.tumor[i, ] %in% LUAD.G2
  temp[6] <- T

  sample.info <- LUAD.batch.info.tumor[i, ]
  sample.info[!temp] <- ""

  unlist(sample.info)
}))
LUAD.raw.intensity.G2 <- LUAD.raw.intensity.tumor[, as.vector(t(LUAD.batch.info.G2))[nchar(t(LUAD.batch.info.G2)) > 0]]

# Extract batch information and intensity for samples in Group3.
LUAD.batch.info.G3 <- t(sapply(1:nrow(LUAD.batch.info.tumor), function(i) {
  temp <- LUAD.batch.info.tumor[i, ] %in% LUAD.G3
  temp[6] <- T

  sample.info <- LUAD.batch.info.tumor[i, ]
  sample.info[!temp] <- ""

  unlist(sample.info)
}))
LUAD.raw.intensity.G3 <- LUAD.raw.intensity.tumor[, as.vector(t(LUAD.batch.info.G3))[nchar(t(LUAD.batch.info.G3)) > 0]]

# Construct a proObj for the Group1 of tumor samples.
LUAD.tumor.G1 <- proObjFromTMT(LUAD.raw.intensity.G1, LUAD.batch.info.G1,
                               IDs = rownames(LUAD.raw.intensity.G1),
                               name = "LUAD.tumor.G1")

# Perform technical & biological variance estimation procedure.
LUAD.tumor.G1 <- estTechVar(LUAD.tumor.G1)
LUAD.tumor.G1 <- estBioVar(LUAD.tumor.G1)

# Construct a proObj for the Group1 of tumor samples.
LUAD.tumor.G2 <- proObjFromTMT(LUAD.raw.intensity.G2, LUAD.batch.info.G2,
                               IDs = rownames(LUAD.raw.intensity.G2),
                               name = "LUAD.tumor.G2")

# Perform technical & biological variance estimation procedure.
LUAD.tumor.G2 <- estTechVar(LUAD.tumor.G2)
LUAD.tumor.G2 <- estBioVar(LUAD.tumor.G2)

# Construct a proObj for the Group1 of tumor samples.
LUAD.tumor.G3 <- proObjFromTMT(LUAD.raw.intensity.G3, LUAD.batch.info.G3,
                               IDs = rownames(LUAD.raw.intensity.G3),
                               name = "LUAD.tumor.G3")

# Perform technical & biological variance estimation procedure.
LUAD.tumor.G3 <- estTechVar(LUAD.tumor.G3)
LUAD.tumor.G3 <- estBioVar(LUAD.tumor.G3)

# Conduct hypothesis testing procedures for identifying DEPs between Group1 & Group2.
DEP.12.G1High.Res <- DiffMean(LUAD.tumor.G1, LUAD.tumor.G2, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")
DEP.12.G2High.Res <- DiffMean(LUAD.tumor.G2, LUAD.tumor.G1, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")

# Conduct hypothesis testing procedures for identifying DEPs between Group1 & Group3.
DEP.13.G1High.Res <- DiffMean(LUAD.tumor.G1, LUAD.tumor.G3, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")
DEP.13.G3High.Res <- DiffMean(LUAD.tumor.G3, LUAD.tumor.G1, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")

# Conduct hypothesis testing procedures for identifying DEPs between Group2 & Group3.
DEP.23.G2High.Res <- DiffMean(LUAD.tumor.G2, LUAD.tumor.G3, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")
DEP.23.G3High.Res <- DiffMean(LUAD.tumor.G3, LUAD.tumor.G2, scaleV.raw.FC = 0,
                              alternative = "upper", p.adjust.method = "BH")

# Identify DEPs under the condition that the BH-corrected P-value is < 0.01.
DEP.G1High <- intersect(rownames(DEP.12.G1High.Res)[which(DEP.12.G1High.Res$padj < 0.01)],
                        rownames(DEP.13.G1High.Res)[which(DEP.13.G1High.Res$padj < 0.01)])
DEP.G2High <- intersect(rownames(DEP.12.G2High.Res)[which(DEP.12.G2High.Res$padj < 0.01)],
                        rownames(DEP.23.G2High.Res)[which(DEP.23.G2High.Res$padj < 0.01)])
DEP.G3High <- intersect(rownames(DEP.13.G3High.Res)[which(DEP.13.G3High.Res$padj < 0.01)],
                        rownames(DEP.23.G3High.Res)[which(DEP.23.G3High.Res$padj < 0.01)])

# Survival analysis -------------------------------------------------------

# Extract overall survival information of LUAD patients.
LUAD.surv.OS <- data.frame(time = LUAD.clinical.info[, "OS_days"],
                           status = LUAD.clinical.info[, "OS_event"],
                           row.names = rownames(LUAD.clinical.info))
LUAD.surv.OS$subgroup <- LUAD.tumor.DVP.3CRes[rownames(LUAD.surv.OS)]

# Create survival curve.
LUAD.surv.fit.Res <- survfit(Surv(time, status) ~ subgroup,
                             data = LUAD.surv.OS)

# Compute p-value by using the log-rank test.
LUAD.surv.diff.Pval <- signif(surv_pvalue(LUAD.surv.fit.Res, data = LUAD.surv.OS)$pval, 3)
Leg <- sapply(1:3, function(i)
  paste0("G", i, " (", length(which(LUAD.surv.OS$subgroup == i)),
         " patients)"))

# Plot Kaplan-Meier survival curve.
LUAD.sur.curve <- ggsurvplot(
  fit = LUAD.surv.fit.Res, data = LUAD.surv.OS,
  pval = LUAD.surv.diff.Pval,
  pval.coord = c(50*0.6, 0.1),
  pval.size = 8, size = 2.2,
  xlab = "Overall survival (day)", ylab = "Survival ratio",
  title = "740 DVPs (BH-adjusted p < 0.01)",
  xlim = c(0, 2000),
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
clus.order <- LUAD.tumor.DVP.CCRes[[3]]$consensusTree$order
bk <- seq(-0.6, 0.6, by = 0.05)

LUAD.clinical.category <- LUAD.clinical.info[, 1:4]
LUAD.clinical.category$`STK11.mutation` <- toupper(LUAD.clinical.info$STK11_mutation)
LUAD.clinical.category$`KRAS.mutation` <- toupper(LUAD.clinical.info$KRAS_mutation)
LUAD.clinical.category$`EGFR.mutation` <- toupper(LUAD.clinical.info$EGFR_mutation)

# Specify the annotations of columns.
anno.col.clinic <- data.frame("STK11 mutation" = as.factor(LUAD.clinical.category$`STK11.mutation`),
                              "KRAS mutation" = as.factor(LUAD.clinical.category$`KRAS.mutation`),
                              "EGFR mutation" = as.factor(LUAD.clinical.category$`EGFR.mutation`),
                              ConClus = as.factor(LUAD.tumor.DVP.3CRes),
                              row.names = names(LUAD.tumor.DVP.3CRes),
                              check.names = F)
anno.col.clinic <- anno.col.clinic[clus.order, ]
anno.col.clinic <- anno.col.clinic[order(anno.col.clinic$ConClus), ]

anno.col.clinic[which(is.na(anno.col.clinic[, 1])), 1] <- "NO"
anno.col.clinic[which(is.na(anno.col.clinic[, 2])), 2] <- "NO"
anno.col.clinic[which(is.na(anno.col.clinic[, 3])), 3] <- "NO"

# Specify annotation_col track colors.
anno.color.clinic <- list(
  "STK11 mutation" = c("YES" = "black", "NO" = "gray"),
  "KRAS mutation" = c("YES" = "black", "NO" = "gray"),
  "EGFR mutation" = c("YES" = "black", "NO" = "gray"),
  ConClus = c("1" = "#FED439FF", "2" = "#709AE1FF", "3" = "#FD7446FF"))

# Plot heatmap of DEPs between 2 subgroups.
LUAD.DEP.heatmap <- pheatmap(
  PostZScrImp.LUAD.tumor[c(DEP.G1High, DEP.G2High, DEP.G2High), rownames(anno.col.clinic)],
  show_rownames = F, show_colnames = F,
  cluster_cols = F, cluster_row = F,
  gaps_col = length(LUAD.G1), gaps_row = length(DEP.G1High),
  gaps_col = c(length(LUAD.G1), length(LUAD.G1)+length(LUAD.G2)),
  gaps_row = c(length(DEP.G1High), length(DEP.G1High)+length(DEP.G2High)),
  main = paste0("G1 (", length(LUAD.G1), " patients), ",
                "G2 (", length(LUAD.G2), " patients) & ",
                "G3 (", length(LUAD.G3), " patients)"),
  annotation_col = anno.col.clinic, annotation_colors = anno.color.clinic,
  col = c(colorRampPalette(colors = c("blue", "white"))(round(length(bk)/2)),
          colorRampPalette(colors = c("white", "red"))(ceiling(length(bk)/2))),
  breaks = bk)


# Enrichment analysis -----------------------------------------------------

# KEGG enrichment analysis of Group1-high DEPs.
DEP.G1High.ID <- bitr(DEP.G1High,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

DEP.G1High.Enrich.GO <- enrichGO(gene = DEP.G1High.ID$ENTREZID,
                                 OrgDb = "org.Hs.eg.db", ont = "ALL",
                                 pvalueCutoff = 0.5, qvalueCutoff = 0.5)

# KEGG enrichment analysis of Group2-high DEPs.
DEP.G2High.ID <- bitr(DEP.G2High,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

DEP.G2High.Enrich.GO <- enrichGO(gene = DEP.G2High.ID$ENTREZID,
                                 OrgDb = "org.Hs.eg.db", ont = "ALL",
                                 pvalueCutoff = 0.5, qvalueCutoff = 0.5)

# KEGG enrichment analysis of Group3-high DEPs.
DEP.G3High.ID <- bitr(DEP.G3High,
                      fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db")

DEP.G3High.Enrich.GO <- enrichGO(gene = DEP.G3High.ID$ENTREZID,
                                 OrgDb = "org.Hs.eg.db", ont = "ALL",
                                 pvalueCutoff = 0.5, qvalueCutoff = 0.5)


# Feature association evaluating ------------------------------------------

# Evaluating the association between protein expression and patient prognosis.

# Fit a Cox proportional hazards regression model for each protein,
# and the two-tailed p-value produced by it serves as the statistical measure of association significance.
PostM.LUAD.tumor.Cox.Res <- exprCox(
  expr = PostM.LUAD.tumor[, rownames(LUAD.clinical.info)],
  time = LUAD.clinical.info$OS_days,
  status = LUAD.clinical.info$OS_event)

# Evaluating the association between protein expression and somatic mutation status of EGFR.

# Conduct Welch’s one-way test or ANOVA for categorical feature,
# and the resulting p-value serves as the statistical measure of protein-feature association.
PostM.LUAD.tumor.EGFR.Res <- exprCatFeature(expr = PostM.LUAD.tumor[, rownames(LUAD.clinical.category)],
                                            feature = LUAD.clinical.category[, "EGFR.mutation"])

# Evaluating the association between protein expression and somatic mutation status of KRAS.
PostM.LUAD.tumor.KRAS.Res <- exprCatFeature(expr = PostM.LUAD.tumor[, rownames(LUAD.clinical.category)],
                                            feature = LUAD.clinical.category[, "KRAS.mutation"])

# Evaluating the association between protein expression and somatic mutation status of STK11.
PostM.LUAD.tumor.STK11.Res <- exprCatFeature(expr = PostM.LUAD.tumor[, rownames(LUAD.clinical.category)],
                                             feature = LUAD.clinical.category[, "STK11.mutation"])

