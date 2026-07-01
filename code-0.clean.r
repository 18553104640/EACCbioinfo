rm(list = ls())
library(GEOquery)
library(org.Hs.eg.db)
library(tidyverse)
library(ggsci)
library(ggrepel)
library(AnnoProbe)

df.count = read.delim2("rawdata/readcount.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.fpkm = read.delim2("rawdata/fpkm_genename.xls", row.names = 1)
df.annot = df.fpkm[21:27]
df.fpkm = df.fpkm[-(21:27)]
df.annot$geneid = mapIds(org.Hs.eg.db, df.annot$Gene.name, "ENTREZID", "SYMBOL",
                         multiVals = "first")
df.annot = df.annot[!is.na(df.annot$geneid),]
df.annot = df.annot[!duplicated(df.annot$geneid),]

df.count = df.count[rownames(df.annot),]
df.fpkm = df.fpkm[rownames(df.annot),]
rownames(df.count) = df.annot$Gene.name
rownames(df.fpkm) = df.annot$Gene.name
df.count = cbind(gene = rownames(df.count), df.count)
df.fpkm = cbind(gene = rownames(df.fpkm), df.fpkm)

df.count = df.count[-which(colnames(df.count)=="Diease6")]
df.fpkm = df.fpkm[-which(colnames(df.fpkm)=="Diease6")]
colnames(df.count) = str_replace(colnames(df.count), "Diease", "Disease")
colnames(df.fpkm) = str_replace(colnames(df.fpkm), "Diease", "Disease")

write.table(df.count, "cleandata/count.xls", quote = F, sep = "\t", row.names = F)
write.table(df.fpkm, "cleandata/fpkm.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.fpkm = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> 
  as.data.frame()
df.fpkm = cbind(gene = rownames(df.fpkm), df.fpkm)
df.fpkm = reshape2::melt(df.fpkm, id.var = "gene")
df.fpkm$group = str_remove_all(df.fpkm$variable, "[0-9]")
ggplot(df.fpkm, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_manual(values = pal_jama()(2)[2:1], guide = guide_none()) +
  theme_bw() + 
  xlab("Sample") + ylab("log2 FPKM") +
  theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1))
ggsave("01.QC/02.FPKM.Box.png", width = 7, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("01.QC/02.FPKM.Box.pdf", width = 7, height = 4, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.fpkm = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> 
  t() |> as.data.frame()
df.fpkm = df.fpkm[colSums(df.fpkm)>0]
res.pca = prcomp(df.fpkm, center = T, scale. = T)
df.pca = predict(res.pca)[,1:2] |> as.data.frame()
df.pca$sample = rownames(df.pca)
df.pca$group = str_remove_all(df.pca$sample, "[0-9]")
df.pca$group = factor(df.pca$group, levels = c("NC","Disease"))

ggplot(df.pca, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = group)) +
  geom_text_repel(aes(label = sample), size = 3) +
  stat_ellipse(aes(group = group, color = group, fill = group), alpha = 0.7) +
  scale_color_jama() +
  scale_fill_jama() +
  theme_bw() + 
  theme(legend.position = "inside", legend.position.inside = c(0.15,0.85), 
        legend.title = element_blank(), aspect.ratio = 1)
ggsave("01.QC/03.PCA.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("01.QC/03.PCA.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
gs = getGEO(filename = "rawdata/GSE102673_series_matrix.txt.gz", getGPL = F)
df.exp = gs |> exprs() |> as.data.frame()
df.pheno = gs |> pData() |> as.data.frame()
df.annot = idmap("GPL21827", type = "pipe")
df.exp = cbind(probe_id = rownames(df.exp), df.exp)
df.exp = merge(df.annot, df.exp, all = T)
df.exp = df.exp[-1]
colnames(df.exp)[1] = "gene"
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
write.table(df.exp, "cleandata/GSE102673.exp.xls", sep = "\t", quote = F, row.names = F)
write.table(df.pheno, "cleandata/GSE102673.pheno.xls", sep = "\t", quote = F, row.names = F)
