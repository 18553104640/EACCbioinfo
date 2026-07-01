rm(list = ls())
library(GSVA)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(corrplot)
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)

gs = readxl::read_excel("../database/PCD.gene.xlsx")
gs = cbind(id = rownames(gs), gs)

gs = reshape2::melt(gs, id.var = "id")
gs = na.omit(gs[2:3])
gs = split(gs$value, gs$variable)
names(gs) = str_to_sentence(names(gs))

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.pheno = data.frame(sample = colnames(df.exp), 
                      group = str_remove_all(colnames(df.exp), "[0-9]"))

gsvaPar = ssgseaParam(as.matrix(df.exp), gs)
res.gsva = gsva(gsvaPar)
res.gsva = t(res.gsva) |> as.data.frame()
res.gsva$sample = rownames(res.gsva)
df.full = merge(df.pheno, res.gsva)
write.table(df.full, "04.PCD/01.PCD.xls", quote = F, sep = "\t", row.names = F)

p.diff = apply(df.full[-c(1:2)], 2, \(x){
  wilcox.test(x~df.full$group)$p.value
})
df.diff = data.frame(PCD = colnames(df.full)[-c(1:2)],
                     p.val = p.diff,
                     p.adj = p.adjust(p.diff, method = "BH"))
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
write.table(df.diff, "04.PCD/02.PCD.Diff.xls", quote = F, sep = "\t", row.names = F)

df.diff = subset(df.diff, p.adj < 0.05)
df.diff$lab = ifelse(df.diff$p.adj > 0.01, "*",
                     ifelse(df.diff$p.adj > 0.001, "**", "***"))
df.plot = reshape2::melt(df.full, id.var = c("sample","group"), 
                         variable.name = "PCD")
df.plot$lab = df.diff$lab[match(df.plot$PCD, df.diff$PCD)]
df.plot$lab = paste0(str_to_sentence(df.plot$PCD), "\n", df.plot$lab)
df.plot = subset(df.plot, PCD %in% df.diff$PCD)
df.plot$PCD = factor(df.plot$PCD, levels = df.diff$PCD)
df.plot = df.plot[order(df.plot$PCD, decreasing = F),]
df.plot$lab = factor(df.plot$lab, levels = unique(df.plot$lab))
df.plot$group = factor(df.plot$group, levels = c("NC","Disease"))

ggplot(df.plot, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group, color = group), alpha = 0.5, linewidth = 1) +
  geom_boxplot(aes(color = group), fill = "white", width = 0.2) +
  scale_fill_jama() + 
  scale_color_jama() +
  xlab("Group") + ylab("ssGSEA Score") +
  facet_wrap(~lab, scales = "free", ncol = 3) + 
  guides(fill = guide_none(), color = guide_none()) +
  theme_bw()
ggsave("04.PCD/03.PCD.Diff.png", width = 7, height = 7, units = "in", bg = "white", dpi = 300)
ggsave("04.PCD/03.PCD.Diff.pdf", width = 7, height = 7, units = "in", bg = "white", dpi = 300)

pcd = df.diff$PCD
gs = gs[pcd]
gs = purrr::reduce(gs,c) |> unique()
write.table(data.frame(gene = gs),
            "04.PCD/04.PCD.Diff.Gene.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.deg = read.delim2("02.DEG/09.Co-DEG.xls")
df.pcd = read.delim2("04.PCD/04.PCD.Diff.Gene.xls")
df.wgcna = read.delim2("03.WGCNA/05.GeneModule.xls")
df.wgcna = subset(df.wgcna, Module %in% c("brown","yellow"))

ggvenn(list(DEG = df.deg$gene, WGCNA = df.wgcna$Gene, PCD = df.pcd$gene), 
       stroke_size = 0.3, set_name_size = 7, text_size = 5, show_percentage = F) + 
  ylim(c(-2,2)) + 
  scale_fill_bmj()
ggsave("05.Gene/01.Venn.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Gene/01.Venn.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

gs = intersect(df.deg$gene, df.wgcna$Gene) |> intersect(df.pcd$gene)
gs = unique(gs)
write.table(data.frame(gene = gs), "05.Gene/02.Intersection.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
library(enrichplot)
df.deg = read.delim2("05.Gene/02.Intersection.xls")
gs = unique(df.deg$gene)
df.annot = AnnotationDbi::select(org.Hs.eg.db, gs, "ENTREZID", "SYMBOL")

res.go = enrichGO(df.annot$ENTREZID, org.Hs.eg.db, ont = "ALL", readable = T)
res.kegg = enrichKEGG(df.annot$ENTREZID, keyType = "ncbi-geneid")
res.kegg = setReadable(res.kegg, org.Hs.eg.db, "ENTREZID")

df.go = res.go@result
write.table(df.go, "05.Gene/03.GO.xls", quote = F, sep = "\t", row.names = F)

dotplot(res.go, showCategory = 15, label_format = 100)
ggsave("05.Gene/04.GO.png", width = 9, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("05.Gene/04.GO.pdf", width = 9, height = 6, units = "in", dpi = 300, bg = "white")

