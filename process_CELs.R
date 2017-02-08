library(affy)
library(annotate)
library(hgu133plus2.db)
library(sva)
library(xlsx)
#library(R.devices)
#library(aroma.affymetrix)
#library("gdata")



setwd('/Users/moshesilverstein/Desktop/cmap')

#csR <- AffymetrixCelSet$byName('/Users/moshesilverstein/Desktop/cmap/cmap_build02.volume1of7', chipType="HT_HG-U133A")

#Data <- ReadAffy()
#eset <- expresso(Data,bgcorrect.method="mas",
#                 normalize.method="qspline",
#                 pmcorrect.method="subtractmm",
#                 summary.method="playerout")

#eset <- expresso(Data, bgcorrect.method="rma",
#                  normalize.method="quantiles",
#                 pmcorrect.method="pmonly",
#                 summary.method="medianpolish")
eset <- read.csv('/Users/moshesilverstein/Desktop/cmap/cmap_merged.csv')

probes <- rownames(exprs(eset))

genes_df <- select(hgu133plus2.db, probes, c("SYMBOL"), keytype="PROBEID")

# dimwrite.csv(exprs(eset), 'ExprMat_mas5.csv')
#dimwrite.csv(exprs(eset), 'ExprMat_rma.csv')
#write.csv(exprs(eset), 'ExprMat_rma.csv')
write.csv(genes_df, 'Probes2Genes.csv')

## Apply ComBat to correct for batch effect
meta_df <- read.xlsx('/Users/moshesilverstein/Desktop/cmap/cmap_instances_02.xlsx', sheetName = "build02")
#meta_df <- read.xls('cmap_instances_02.xls', sheetName = "Sheet1")
#meta_df <- read.table('/Users/moshesilverstein/Desktop/cmap/cmap_instances_02.xls')

# order meta_df by column order in eset
meta_df <- meta_df[colnames(eset),]
#modcombat <- model.matrix(~1, data=meta_df)
modcombat <- model.matrix(~1, data=eset)

combat_edata <- ComBat(dat=exprs(eset), batch=meta_df$batch_id, 
                       mod=modcombat, 
                       par.prior=TRUE, prior.plots=FALSE)
                       
#combat_edata <- ComBat(dat=eset, batch=meta_df$batch_id, 
#                       mod=modcombat, 
#                       par.prior=TRUE, prior.plots=FALSE)

write.csv(combat_edata, 'ExprMat_rma_combat.csv')

## Average probes
genes_df <- genes_df[!is.na(genes_df$SYMBOL),]
mat <- merge(as.data.frame(exprs(eset)), genes_df, 
             by.x = 'row.names',
             by.y = 'PROBEID',
             all.x = F, all.y = F)

mat['Row.names'] <- NULL
symbols <- mat$SYMBOL
mat['SYMBOL'] <- NULL
mat <- aggregate(mat, by = list(symbols), mean)
rownames(mat) <- mat$Group.1
mat['Group.1'] <- NULL

write.csv(mat, 'ExprMat_rma_avg_probes.csv')

