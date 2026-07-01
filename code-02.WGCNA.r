rm(list = ls())
library(tidyverse)
library(ggsci)
library(WGCNA)

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> t() |> as.data.frame()
df.exp = df.exp[colSums(df.exp)>0]
df.pheno = data.frame(sample = rownames(df.exp), 
                      group = str_remove_all(rownames(df.exp), "[0-9]"))

goodSamplesGenes(df.exp, verbose = 3)
collectGarbage()

par(mfrow = c(1,1))
sampleTree = hclust(dist(df.exp), method = "average")

png("03.WGCNA/01.OutlierDetection.png", width = 7, height = 5, res = 1200, units = "in")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

pdf("03.WGCNA/01.OutlierDetection.pdf", width = 7, height = 5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

powers = c(1:20)
collectGarbage()
sft = pickSoftThreshold(df.exp, powerVector = powers, verbose = 5, blockSize = ncol(df.exp))

png("03.WGCNA/02.SoftThred.png", width = 10, height = 6, units = "in", res = 1200)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

pdf("03.WGCNA/02.SoftThred.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate

cor = WGCNA::cor
net = blockwiseModules(
  df.exp, power = sft$powerEstimate, maxBlockSize = ncol(df.exp), 
  minModuleSize = 100, mergeCutHeight = 0.15, deepSplit = 4, corType = "bicor",
  pamRespectsDendro = FALSE, reassignThreshold = 0, 
  robustY = FALSE, numericLabels = TRUE, verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)

png("03.WGCNA/03.Modules.png", width = 12, height = 6, res = 300, units = "in")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F)
dev.off()
pdf("03.WGCNA/03.Modules.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = F)
dev.off()

moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(df.exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
design = df.pheno
design$group = factor(design$group, levels = c("NC","Disease"))

moduleTraitCor = cor(MEs, as.numeric(design$group))
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(df.exp)) %>% signif(3)

sizeGrWindow(15,10)
textMatrix = paste(signif(moduleTraitCor, 3), "\n(",
                   signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

png("03.WGCNA/04.ModuleTraitCorrelation.png", width = 5, height = 9, units = "in", res = 300)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 'Disease',
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs) %>% str_remove("ME"),
               colorLabels = FALSE,
               keepLegendSpace = F,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               xLabelsAngle = 0, xLabelsAdj = 0.5,
               setStdMargins = T,
               cex.text = 0.7,
               zlim = c(-1,1))
dev.off()

pdf("03.WGCNA/04.ModuleTraitCorrelation.pdf", width = 5, height = 9)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 'Disease',
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs) %>% str_remove("ME"),
               colorLabels = FALSE,
               keepLegendSpace = F,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               xLabelsAngle = 0, xLabelsAdj = 0.5,
               setStdMargins = T,
               cex.text = 0.7,
               zlim = c(-1,1))
dev.off()

genes.color = data.frame(names(net$colors),labels2colors(net$colors))
colnames(genes.color) = c("Gene", "Module")
table(genes.color$Module)
write.table(genes.color, "03.WGCNA/05.GeneModule.xls", sep = "\t", quote = F, col.names = T, row.names = F)





