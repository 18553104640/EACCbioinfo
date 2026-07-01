rm(list = ls())
library(tidyverse)
library(ggsci)
library(ggpubr)
library(IOBR)
library(GSVA)
library(ggrepel)

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.cb = deconvo_cibersort(df.exp, arrays = F, perm = 1)
df.cb = df.cb[1:23]
write.table(df.cb, "interm/cibersort.tsv", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.cb = read.delim2("interm/cibersort.tsv")
df.cb = reshape2::melt(df.cb, id.var = "ID", variable.name = "cell", value.name = "value")
df.cb$value = as.numeric(df.cb$value)
df.cb$cell = str_remove(df.cb$cell, "_CIBERSORT")
df.cb$group = str_remove_all(df.cb$ID, "[0-9]")
df.cb$group = factor(df.cb$group, levels = c("NC", "Disease"))

cell.cols = palettes(category = "random", 22, show_col = F, show_message = F)[1:22]
ggplot(df.cb, aes(x = ID, y = value, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~group, scale = "free") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = cell.cols, guide = guide_legend(ncol = 1)) +
  xlab("Sample") + ylab("Proportion") +
  theme_bw() + 
  theme(axis.text.x.bottom = element_text(hjust = 1, angle = 30),
        strip.background = element_blank())
ggsave("07.Immune/01.CIBERSORT.png", width = 8, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("07.Immune/01.CIBERSORT.pdf", width = 8, height = 7, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.cb = read.delim2("interm/cibersort.tsv", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.cb = df.cb[colSums(df.cb)>0]
res.pca = prcomp(df.cb, center = T, scale. = T)
df.pca = predict(res.pca) |> as.data.frame()
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
ggsave("07.Immune/02.Immune.PCA.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("07.Immune/02.Immune.PCA.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.cb = read.delim2("interm/cibersort.tsv", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
grp = str_remove_all(rownames(df.cb), "[0-9]")

p.diff = apply(df.cb, 2, \(x){
  wilcox.test(x~grp)$p.value
})
df.diff = data.frame(cell = colnames(df.cb), 
                     p.val = p.diff,
                     p.adj = p.adjust(p.diff, method = "BH"))
write.table(df.diff, '07.Immune/03.Immune.Diff.xls', quote = F, sep = "\t", row.names = F)

df.diff = subset(df.diff, p.adj < 0.05)
df.cb = df.cb[df.diff$cell]
df.cb = cbind(sample = rownames(df.cb), df.cb)
df.cb = reshape2::melt(df.cb, id.var = 'sample', variable.name = "cell")
df.cb$cell = str_remove(df.cb$cell, "_CIBERSORT")
df.diff$cell = str_remove(df.diff$cell, "_CIBERSORT")
df.diff$lab = ifelse(df.diff$p.adj < 0.01, "*", 
                     ifelse(df.diff$p.adj > 0.001, "**", "***"))
df.cb$group = str_remove_all(df.cb$sample, "[0-9]")
df.cb$group = factor(df.cb$group, levels = c("NC","Disease"))

ggplot(df.cb, aes(x = cell, y = value)) + 
  geom_boxplot(aes(color = group)) +
  geom_text(aes(y = 1, label = lab), data = df.diff) +
  scale_color_jama() +
  xlab(NULL) + ylab("Proportion") +
  theme_classic2() + 
  theme(legend.position = "top", legend.title = element_blank())
ggsave("07.Immune/04.Immune.Diff.png", width = 6, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("07.Immune/04.Immune.Diff.pdf", width = 6, height = 6, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.cb = read.delim2("interm/cibersort.tsv", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
colnames(df.cb) = str_remove(colnames(df.cb), "_CIBERSORT")
df.cb = df.cb[str_starts(colnames(df.cb), "Macro")]

df.cor = cor(df.cb, method = "spearman")
df.p = cor.mtest(df.cb, method = "spearman")$p

df.cor.out = cbind(cell1 = rownames(df.cor), as.data.frame(df.cor))
df.p.out = cbind(cell1 = rownames(df.p), as.data.frame(df.p))
df.cor.out = reshape2::melt(df.cor.out, id.var = "cell1", variable.name = "cell2", value.name = "cor")
df.p.out = reshape2::melt(df.p.out, id.var = "cell1", variable.name = "cell2", value.name = "p")
df.out = merge(df.cor.out, df.p.out)
write.table(df.out, "07.Immune/05.Immune.Diff.Cor.xls", quote = F, sep = "\t", row.names = F)

png("07.Immune/06.Immune.Diff.Cor.png", width = 5, height = 5, units = "in", bg = "white", res = 300)
corrplot(df.cor, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = df.p, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
dev.off()
pdf("07.Immune/06.Immune.Diff.Cor.pdf", width = 5, height = 5)
corrplot(df.cor, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = df.p, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
dev.off()


rm(list = ls())
df.cb = read.delim2("interm/cibersort.tsv", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
colnames(df.cb) = str_remove(colnames(df.cb), "_CIBERSORT")
df.cb = df.cb[str_starts(colnames(df.cb), "Macro")]

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.exp = t(df.exp) |> as.data.frame()
df.sel = read.delim2("06.Selection/06.Selected.xls")
df.exp = df.exp[rownames(df.cb), df.sel$gene]

res.cor = linkET::correlate(df.cb, df.exp, method = "spearman")
df.cor = res.cor$r
df.p = res.cor$p

df.cor.out = cbind(cell = rownames(df.cor), as.data.frame(df.cor))
df.p.out = cbind(cell = rownames(df.p), as.data.frame(df.p))
df.cor.out = reshape2::melt(df.cor.out, id.var = "cell", variable.name = "gene", value.name = "cor")
df.p.out = reshape2::melt(df.p.out, id.var = "cell", variable.name = "gene", value.name = "p")
df.out = merge(df.cor.out, df.p.out)
write.table(df.out, "07.Immune/07.Immune.Gene.Cor.xls", quote = F, sep = "\t", row.names = F)

ggplot(df.out, aes(x = gene, y = cell)) + 
  geom_point(aes(color = cor, size = abs(cor))) + 
  geom_text(aes(label = signif(cor, 2)), size = 3) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
  xlab(NULL) + ylab(NULL) +
  scale_size_continuous(range = c(0,9), guide = guide_none(), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1))
ggsave("07.Immune/08.Immune.Gene.Cor.png", width = 6, height = 3, units = "in", dpi = 300, bg = "white")
ggsave("07.Immune/08.Immune.Gene.Cor.pdf", width = 6, height = 3, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x){log2(as.numeric(x)+1)}) |> as.data.frame()

gs = read.delim2("../database/TIP.txt")
gs$set = paste0(gs$Steps,gs$Direction)
gs = split(gs$GeneSymbol, gs$set)

in.gsva = gsvaParam(as.matrix(df.exp), gs)
res.gsva = gsva(in.gsva) |> as.data.frame()
res.gsva = t(res.gsva) |> as.data.frame()

res.gsva$step_1 = res.gsva$`1positive`
res.gsva$step_2 = res.gsva$`2positive` - res.gsva$`2negative`
res.gsva$step_3 = res.gsva$`3positive` - res.gsva$`3negative`
res.gsva$step_4 = res.gsva$`4positive`
res.gsva$step_5 = res.gsva$`5positive` - res.gsva$`5negative`
res.gsva$step_6 = res.gsva$`6positive` - res.gsva$`6negative`
res.gsva$step_7 = res.gsva$`7positive` - res.gsva$`7negative`
res.gsva = res.gsva[str_starts(colnames(res.gsva), "step")]
res.gsave = res.gsva[str_starts(rownames(res.gsva), "Disease"),]

df.exp = t(df.exp) |> as.data.frame()
df.exp = df.exp[rownames(res.gsva),]
df.sel = read.delim2("06.Selection/06.Selected.xls")
df.exp = df.exp[df.sel$gene]

res.cor = linkET::correlate(res.gsva, df.exp, method = "spearman")
df.cor = res.cor$r
df.p = res.cor$p

df.cor.out = cbind(step = rownames(df.cor), as.data.frame(df.cor))
df.p.out = cbind(step = rownames(df.p), as.data.frame(df.p))
df.cor.out = reshape2::melt(df.cor.out, id.var = "step", variable.name = "gene", value.name = "cor")
df.p.out = reshape2::melt(df.p.out, id.var = "step", variable.name = "gene", value.name = "p")
df.out = merge(df.cor.out, df.p.out)

name.dict = c(
  step_1 = "Step 1 Release of cancer cell antigens",
  step_2 = "Step 2 Cancer antigen presentation",
  step_3 = "Step 3 Priming and activation",
  step_4 = "Step 4 Trafficking of immune cells to tumors",
  step_5 = "Step 5 Infiltration of immune cells into tumors",
  step_6 = "Step 6 Recognition of cancer cells by T cells",
  step_7 = "Step 7 Killing of cancer cells"
)

df.out$step = name.dict[df.out$step]
write.table(df.out, "07.Immune/09.TIP.Cor.xls", quote = F, sep = "\t", row.names = F)

df.out$step = factor(df.out$step, levels = rev(name.dict))
ggplot(df.out, aes(x = gene, y = step)) + 
  geom_point(aes(color = cor, size = abs(cor))) + 
  geom_text(aes(label = signif(cor, 2)), size = 3) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
  xlab(NULL) + ylab(NULL) +
  scale_size_continuous(range = c(0,9), guide = guide_none(), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1))
ggsave("07.Immune/10.TIP.Cor.png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("07.Immune/10.TIP.Cor.pdf", width = 7, height = 5, units = "in", dpi = 300, bg = "white")



