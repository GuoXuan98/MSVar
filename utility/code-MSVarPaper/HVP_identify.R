# MSVar identifies HVPs across tumor samples for each data set.

library(MSVar)
library(openxlsx)

# BRCA, tumor -------------------------------------------------------------

# Read data.
BRCA.raw.intensity <- read.xlsx("data/BRCA_122Tumor.xlsx",
                                rowNames = T, check.names = F, sheet = 1)
BRCA.batch.info <- read.xlsx("data/BRCA_122Tumor.xlsx",
                             rowNames = F, check.names = F, sheet = 2)

# Construct a proObj for the group of tumor samples.
BRCA.tumor <- proObjFromTMT(BRCA.raw.intensity, BRCA.batch.info,
                            IDs = rownames(raw.intensity),
                            name = "BRCA.tumor")

# Perform technical variance estimation procedure.
BRCA.tumor <- estTechVar(BRCA.tumor)

# Perform biological variance estimation procedure.
BRCA.tumor <- estBioVar(BRCA.tumor)

# Derive posterior M-value.
BRCA.tumor <- PostM(BRCA.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across BRCA tumor samples.
BRCA.HVP <- VarTest(BRCA.tumor)

# CC, tumor -------------------------------------------------------------

# Read data.
CC.raw.intensity <- read.xlsx("data/CC_136Tumor.xlsx",
                              rowNames = T, check.names = F, sheet = 1)
CC.batch.info <- read.xlsx("data/CC_136Tumor.xlsx",
                           rowNames = F, check.names = F, sheet = 2)

# Construct a proObj for the group of tumor samples.
CC.tumor <- proObjFromTMT(CC.raw.intensity, CC.batch.info,
                          IDs = rownames(raw.intensity),
                          name = "CC.tumor")

# Perform technical variance estimation procedure.
CC.tumor <- estTechVar(CC.tumor)

# Perform biological variance estimation procedure.
CC.tumor <- estBioVar(CC.tumor)

# Derive posterior M-value.
CC.tumor <- PostM(CC.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across CC tumor samples.
CC.HVP <- VarTest(CC.tumor)

# GBM, tumor -------------------------------------------------------------

# Read data.
GBM.raw.intensity <- read.xlsx("data/GBM_99Tumor.xlsx",
                               rowNames = T, check.names = F, sheet = 1)
GBM.batch.info <- read.xlsx("data/GBM_99Tumor.xlsx",
                            rowNames = F, check.names = F, sheet = 2)

# Construct a proObj for the group of tumor samples.
GBM.tumor <- proObjFromTMT(GBM.raw.intensity, GBM.batch.info,
                           IDs = rownames(raw.intensity),
                           name = "GBM.tumor")

# Perform technical variance estimation procedure.
GBM.tumor <- estTechVar(GBM.tumor)

# Perform biological variance estimation procedure.
GBM.tumor <- estBioVar(GBM.tumor)

# Derive posterior M-value.
GBM.tumor <- PostM(GBM.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across GBM tumor samples.
GBM.HVP <- VarTest(GBM.tumor)

# HGSOC, tumor -------------------------------------------------------------

# Read data.
HGSOC.raw.intensity <- read.xlsx("data/HGSOC_82Tumor.xlsx",
                                 rowNames = T, check.names = F, sheet = 1)
HGSOC.batch.info <- read.xlsx("data/HGSOC_82Tumor.xlsx",
                              rowNames = F, check.names = F, sheet = 2)

# Construct a proObj for the group of tumor samples.
HGSOC.tumor <- proObjFromTMT(HGSOC.raw.intensity, HGSOC.batch.info,
                             IDs = rownames(raw.intensity),
                             name = "HGSOC.tumor")

# Perform technical variance estimation procedure.
HGSOC.tumor <- estTechVar(HGSOC.tumor)

# Perform biological variance estimation procedure.
HGSOC.tumor <- estBioVar(HGSOC.tumor)

# Derive posterior M-value.
HGSOC.tumor <- PostM(HGSOC.tumor)

# Conduct hypothesis testing procedures for identifying HVPs across HGSOC tumor samples.
HGSOC.HVP <- VarTest(HGSOC.tumor)
