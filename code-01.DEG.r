rm(list = ls())
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(limma)

df.count = read.delim2("cleandata/count.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
grps = str_remove_all(colnames(df.count), "[0-9]")
df.grp = data.frame(sample = colnames(df.count), Group = grps, 
                    row.names = colnames(df.count))
df.grp$Group = factor(df.grp$Group, levels = c("NC","Disease"))

dds = DESeqDataSetFromMatrix(df.count, df.grp, design = ~Group)
dds = DESeq(dds)
res.dds = results(dds) %>% as.data.frame
res.dds = na.omit(res.dds)
res.dds = res.dds[c(2,5,6)]
colnames(res.dds) = c("logFC","P.Value","adj.P.Val")
res.dds = cbind(Gene = rownames(res.dds), res.dds)

res.dds$Direction = ifelse(res.dds$logFC > 1, "Up", "No")
res.dds$Direction = ifelse(res.dds$logFC < -1, "Down", res.dds$Direction)
res.dds$Direction = ifelse(res.dds$adj.P.Val < 0.05, res.dds$Direction, "No")
res.dds.out = subset(res.dds, Direction != "No")
write.table(res.dds.out, "02.DEG/01.RNAseq.DEG.xls", quote = F, sep = "\t", row.names = F)

ggplot(res.dds, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  scale_color_manual(values=c("steelblue", "grey30", "salmon"),
                     labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  geom_vline(xintercept = c(-1,1), col="grey60", linetype = 2) +
  scale_x_continuous(limits = c(-28,28)) +
  ggtitle("Disease vs NC") +
  ylab("-log10(p.adj)") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
ggsave("02.DEG/02.RNAseq.Volcano.png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/02.RNAseq.Volcano.pdf", width = 7, height = 5, units = "in", dpi = 300, bg = "white")

res.dds = res.dds[order(res.dds$Direction, res.dds$adj.P.Val, decreasing = F),]
res.dds.sub = rbind(subset(res.dds, Direction == "Up")[1:10,],
                    subset(res.dds, Direction == "Down")[1:10,])
df.fpkm = read.delim2("cleandata/fpkm.xls", row.names = 1)
df.sub = df.fpkm[rownames(res.dds.sub),]
col.mat = data.frame(Sample = df.grp$sample, Group = df.grp$Group)
col.mat = col.mat[order(col.mat$Group),]
rownames(col.mat) = col.mat$Sample
row.mat = data.frame(Direction = factor(res.dds.sub$Direction, levels = c("Up", "Down")))
rownames(row.mat) = rownames(df.sub)
df.sub = df.sub[rownames(col.mat)]
annotation.cols = list(Group = c(`NC` = "#374E55", `Disease` = "#DF8F44"),
                       Direction = c(Up = "red", Down = "blue"))
df.sub = apply(df.sub, c(1,2), \(x) log2(as.numeric(x)+1))
p = pheatmap(df.sub, scale = "row", border_color = "white",
             cluster_rows = T, cluster_cols = F, show_colnames = F, 
             annotation_colors = annotation.cols,
             annotation_names_col = F, annotation_names_row = F,
             annotation_col = col.mat[-1], annotation_row = row.mat)
png("02.DEG/03.RNAseq.Heatmap.png", width = 7, height = 5, units = "in", bg = "white", res = 300)
print(p)
dev.off()
pdf("02.DEG/03.RNAseq.Heatmap.pdf", width = 7, height = 5)
print(p)
dev.off()


rm(list = ls())
df.exp = read.delim2("cleandata/GSE102673.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE102673.pheno.xls")
df.pheno$group = ifelse(df.pheno$tissue.ch1 == "normal skin", "NC", "Disease")
df.exp = df.exp[df.pheno$geo_accession]
design.mat = cbind(Control = ifelse(df.pheno$group == "NC", 1, 0), 
                   Patient = ifelse(df.pheno$group == "NC", 0, 1))
contrast.mat = makeContrasts(contrasts="Patient-Control", levels=design.mat)

fit = lmFit(df.exp, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
fit = fit[c(1,4,5)]
fit = cbind(Gene = rownames(fit), fit)
fit$Direction = ifelse(fit$logFC > 0, "Up", "No")
fit$Direction = ifelse(fit$logFC < 0, "Down", fit$Direction)
fit$Direction = ifelse(fit$P.Value < 0.05, fit$Direction, "No")
fit$Direction = factor(fit$Direction, levels = c("Down", "No", "Up"))
fit = fit[order(fit$adj.P.Val, decreasing = F), ]
fit.out = subset(fit, Direction != "No")
write.table(fit.out, "02.DEG/04.GSE102673.DEG.xls", quote = F, sep = "\t", row.names = F)

ggplot(fit, aes(x = logFC, y = -log10(P.Value), color = Direction)) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  scale_color_manual(values=c("steelblue", "grey30", "salmon"), labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  xlim(-10,10) +
  ggtitle("Disease vs NC") +
  ylab("-log10(pvalue)") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
ggsave("02.DEG/05.GSE102673.Volcano.png", width = 7, height = 5, dpi = 300, bg = "white")
ggsave("02.DEG/05.GSE102673.Volcano.pdf", width = 7, height = 5, dpi = 300, bg = "white")

fit = fit[order(fit$Direction, fit$adj.P.Val, decreasing = F),]
fit = na.omit(fit)
fit.sub = rbind(subset(fit, Direction == "Up")[1:10,],
                subset(fit, Direction == "Down")[1:10,])
df.sub = df.exp[rownames(fit.sub),]
rownames(df.sub) = fit.sub$Gene
col.mat = data.frame(Sample = df.pheno$geo_accession, Group = df.pheno$group)
col.mat = col.mat[order(col.mat$Group),]
rownames(col.mat) = col.mat$Sample
row.mat = data.frame(Direction = factor(fit.sub$Direction, levels = c("Up", "Down")))
rownames(row.mat) = rownames(df.sub)
df.sub = df.sub[rownames(col.mat)]
annotation.cols = list(Group = c(`NC` = "#374E55", `Disease` = "#DF8F44"),
                       Direction = c(Up = "red", Down = "blue"))

p = pheatmap(df.sub, scale = "row", border_color = "white",
             cluster_rows = T, cluster_cols = F, show_colnames = F, 
             annotation_colors = annotation.cols,
             annotation_names_col = F, annotation_names_row = F,
             annotation_col = col.mat[-1], annotation_row = row.mat)
png("02.DEG/06.GSE102673.Heatmap.png", width = 5, height = 5, units = "in", bg = "white", res = 300)
print(p)
dev.off()
pdf("02.DEG/06.GSE102673.Heatmap.pdf", width = 5, height = 5)
print(p)
dev.off()


rm(list = ls())
df.deg.1 = read.delim2("02.DEG/01.RNAseq.DEG.xls")
df.deg.2 = read.delim2("02.DEG/04.GSE102673.DEG.xls")
deg.1.up = subset(df.deg.1, Direction == "Up")$Gene
deg.2.up = subset(df.deg.2, Direction == "Up")$Gene
deg.1.down = subset(df.deg.1, Direction == "Down")$Gene
deg.2.down = subset(df.deg.2, Direction == "Down")$Gene

ggvenn(list(RNAseq = deg.1.up, GSE102673 = deg.2.up), 
       stroke_size = 0.3, set_name_size = 4, text_size = 4, show_percentage = F) + 
  ylim(c(-2,1.5)) + 
  scale_fill_bmj() + 
  ggtitle("Upregulated Genes in Disease") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("02.DEG/07.Co-Up.DEG.png", width = 5, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/07.Co-Up.DEG.pdf", width = 5, height = 4, units = "in", dpi = 300, bg = "white")

ggvenn(list(RNAseq = deg.1.down, GSE102673 = deg.2.down), 
       stroke_size = 0.3, set_name_size = 4, text_size = 4, show_percentage = F) + 
  ylim(c(-2,1.5)) + 
  scale_fill_bmj() + 
  ggtitle("Downregulated Genes in Disease") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("02.DEG/08.Co-Down.DEG.png", width = 5, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/08.Co-Down.DEG.pdf", width = 5, height = 4, units = "in", dpi = 300, bg = "white")

df.up = data.frame(gene = intersect(deg.1.up, deg.2.up), direction = "Up")
df.down = data.frame(gene = intersect(deg.1.down, deg.2.down), direction = "Down")
df.out = rbind(df.up, df.down)
write.table(df.out, "02.DEG/09.Co-DEG.xls", quote = F, sep = "\t", row.names = F)

rm(list = ls())
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

df.deg = read.delim2("02.DEG/09.Co-DEG.xls")
gs = unique(df.deg$gene)
df.annot = AnnotationDbi::select(org.Hs.eg.db, gs, "ENTREZID", "SYMBOL")

res.go = enrichGO(df.annot$ENTREZID, org.Hs.eg.db, ont = "ALL", readable = T)
res.kegg = enrichKEGG(df.annot$ENTREZID, keyType = "ncbi-geneid")
res.kegg = setReadable(res.kegg, org.Hs.eg.db, "ENTREZID")

df.go = res.go@result
write.table(df.go, "02.DEG/10.DEG.GO.xls", quote = F, sep = "\t", row.names = F)

dotplot(res.go, showCategory = 15, label_format = 100)
ggsave("02.DEG/11.DEG.GO.png", width = 9, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/11.DEG.GO.pdf", width = 9, height = 6, units = "in", dpi = 300, bg = "white")

df.kegg = res.kegg@result
write.table(df.kegg, "02.DEG/12.DEG.KEGG.xls", quote = F, sep = "\t", row.names = F)

dotplot(res.kegg, showCategory = 15, label_format = 100)
ggsave("02.DEG/13.DEG.KEGG.png", width = 9, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.DEG/13.DEG.KEGG.pdf", width = 9, height = 4, units = "in", dpi = 300, bg = "white")



