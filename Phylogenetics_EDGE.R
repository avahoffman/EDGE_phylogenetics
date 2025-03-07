SGS <- read.csv("PD_data/SGS/SGS_cover.csv", header = T)
KNZ <- read.csv("PD_data/KNZ/KNZ_cover.csv", header = T)
HPG <- read.csv("PD_data/HPG/HPG_cover.csv", header = T)
HYS <- read.csv("PD_data/HYS/HYS_cover.csv", header = T)
SBL <- read.csv("PD_data/SBL/SBL_cover.csv", header = T)
SBK <- read.csv("PD_data/SBK/SBK_cover.csv", header = T)



library(ape)
library(ggtree)
library(picante)
## look into UNIFRAC?


## SGS 
###########################################################################################

# generate tree from Cipres Data
community_SGS <- read.nexus("PD_data/SGS/gene_data/Cipres_Data/RAxML_bestTree.result")
community_SGS$edge.length <- round(community_SGS$edge.length, digits = 3)
ggtree(community_SGS) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/SGS/figures/Distances_toscale.jpg")

# remove rownames
SGS_noplotnames <- SGS[,-1]

# calculate Faith's PD
tree_SGS <- prune.sample(SGS_noplotnames,community_SGS)
phylo.div_SGS <- pd(SGS_noplotnames,community_SGS, include.root = T)

d <- cbind(SGS[,1],phylo.div_SGS)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/SGS_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_SGS <- ses.pd(SGS_noplotnames,community_SGS,null.model=c("independentswap"))

e <- cbind(SGS[,1],ses.pd_SGS)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/SGS_SESPD.csv")


## KNZ
###########################################################################################

# generate tree from Cipres Data
community_KNZ <- read.nexus("PD_data/KNZ/gene_data/Cipres_Data/RAxML_bestTree.result")
community_KNZ$edge.length <- round(community_KNZ$edge.length, digits = 3)
ggtree(community_KNZ) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/KNZ/figures/Distances_toscale.jpg")

# remove rownames
KNZ_noplotnames <- KNZ[,-1]

# calculate Faith's PD
tree_KNZ <- prune.sample(KNZ_noplotnames,community_KNZ)
phylo.div_KNZ <- pd(KNZ_noplotnames,community_KNZ, include.root = T)

d <- cbind(KNZ[,1],phylo.div_KNZ)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/KNZ_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_KNZ <- ses.pd(KNZ_noplotnames,community_KNZ,null.model=c("independentswap"))

e <- cbind(KNZ[,1],ses.pd_KNZ)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/KNZ_SESPD.csv")


## HPG
###########################################################################################

# generate tree from Cipres Data
community_HPG <- read.nexus("PD_data/HPG/gene_data/Cipres_Data/RAxML_bestTree.result")
community_HPG$edge.length <- round(community_HPG$edge.length, digits = 3)
ggtree(community_HPG) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/HPG/figures/Distances_toscale.jpg")

# remove rownames
HPG_noplotnames <- HPG[,-1]

# calculate Faith's PD
tree_HPG <- prune.sample(HPG_noplotnames,community_HPG)
phylo.div_HPG <- pd(HPG_noplotnames,community_HPG, include.root = T)

d <- cbind(HPG[,1],phylo.div_HPG)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/HPG_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_HPG <- ses.pd(HPG_noplotnames,community_HPG,null.model=c("independentswap"))

e <- cbind(HPG[,1],ses.pd_HPG)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/HPG_SESPD.csv")


## HYS
###########################################################################################

# generate tree from Cipres Data
community_HYS <- read.nexus("PD_data/HYS/gene_data/Cipres_Data/RAxML_bestTree.result")
community_HYS$edge.length <- round(community_HYS$edge.length, digits = 3)
ggtree(community_HYS) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/HYS/figures/Distances_toscale.jpg")

# remove rownames
HYS_noplotnames <- HYS[,-1]

# calculate Faith's PD
tree_HYS <- prune.sample(HYS_noplotnames,community_HYS)
phylo.div_HYS <- pd(HYS_noplotnames,community_HYS, include.root = T)

d <- cbind(HYS[,1],phylo.div_HYS)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/HYS_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_HYS <- ses.pd(HYS_noplotnames,community_HYS,null.model=c("independentswap"))

e <- cbind(HYS[,1],ses.pd_HYS)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/HYS_SESPD.csv")


## SBL
###########################################################################################

# generate tree from Cipres Data
community_SBL <- read.nexus("PD_data/SBL/gene_data/Cipres_Data/RAxML_bestTree.result")
community_SBL$edge.length <- round(community_SBL$edge.length, digits = 3)
ggtree(community_SBL) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/SBL/figures/Distances_toscale.jpg")

# remove rownames
SBL_noplotnames <- SBL[,-1]

# calculate Faith's PD
tree_SBL <- prune.sample(SBL_noplotnames,community_SBL)
phylo.div_SBL <- pd(SBL_noplotnames,community_SBL, include.root = T)

d <- cbind(SBL[,1],phylo.div_SBL)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/SBL_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_SBL <- ses.pd(SBL_noplotnames,community_SBL,null.model=c("independentswap"))

e <- cbind(SBL[,1],ses.pd_SBL)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/SBL_SESPD.csv")


## SBK
###########################################################################################

# generate tree from Cipres Data
community_SBK <- read.nexus("PD_data/SBK/gene_data/Cipres_Data/RAxML_bestTree.result")
community_SBK$edge.length <- round(community_SBK$edge.length, digits = 3)
ggtree(community_SBK) + geom_tiplab() + xlim(NA, 1)
ggsave(width=15,height=8, filename="PD_data/SBK/figures/Distances_toscale.jpg")

# remove rownames
SBK_noplotnames <- SBK[,-1]

# calculate Faith's PD
tree_SBK <- prune.sample(SBK_noplotnames,community_SBK)
phylo.div_SBK <- pd(SBK_noplotnames,community_SBK, include.root = T)

d <- cbind(SBK[,1],phylo.div_SBK)
colnames(d)[1] <- 'plot'
d$pd_per_spp <- d$PD / d$SR

write.csv(d, "PD_final/SBK_PD.csv")

# calculate SESFaith's PD (see Mason et al 2013)
# use default iteration num (1000) and runs (999)
ses.pd_SBK <- ses.pd(SBK_noplotnames,community_SBK,null.model=c("independentswap"))

e <- cbind(SBK[,1],ses.pd_SBK)
colnames(e)[1] <- 'plot'

write.csv(e, "PD_final/SBK_SESPD.csv")