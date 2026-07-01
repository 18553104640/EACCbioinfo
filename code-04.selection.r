rm(list = ls())
library(tidyverse)
library(glmnet)
library(Boruta)
library(caret)
library(ggsci)
library(ggpubr)
library(ggvenn)
library(corrplot)
library(pROC)
library(IOBR)
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
gs = read.delim2("05.Gene/02.Intersection.xls")
gs = gs$gene
df.pheno = data.frame(sample = colnames(df.exp), 
                      group = str_remove_all(colnames(df.exp), "[0-9]"))

df.exp = df.exp[gs, ]
df.exp = t(df.exp) |> as.data.frame()

set.seed(2025)
res.lasso = cv.glmnet(as.matrix(df.exp), df.pheno$group, 
                      family = "binomial", nfolds = 10, alpha = 1)
plot(res.lasso)
ggsave("06.Selection/01.Lasso.CV.png", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("06.Selection/01.Lasso.CV.pdf", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("06.Selection/02.Lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("06.Selection/02.Lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")

df.coef = coef(res.lasso, s = res.lasso$lambda.min)
df.coef = cbind(gene = rownames(df.coef), coefficient = df.coef[,1]) |> as.data.frame()
df.coef$coefficient = as.numeric(df.coef$coefficient)
df.coef = subset(df.coef, coefficient != 0)
df.coef = df.coef[-1,]
df.coef = df.coef[order(df.coef$coefficient, decreasing = F),]
df.coef$gene = factor(df.coef$gene, levels = df.coef$gene)

ggplot(df.coef, aes(y = gene, x = coefficient)) + 
  geom_bar(aes(fill = gene), stat = "identity") + 
  geom_text(aes(label = signif(coefficient,3)), hjust = 0, data = subset(df.coef, coefficient>0)) +
  geom_text(aes(label = signif(coefficient,3)), hjust = 1, data = subset(df.coef, coefficient<0)) +
  scale_fill_bmj() +
  xlim(c(-3.5,2)) + 
  ylab("Gene") + xlab("Coefficient") +
  guides(fill = guide_none()) +
  theme_classic2()
ggsave("06.Selection/03.Lasso.Coefficient.png", width = 6, height = 3, dpi = 300, units = "in", bg = "white")
ggsave("06.Selection/03.Lasso.Coefficient.pdf", width = 6, height = 3, dpi = 300, units = "in", bg = "white")

gs.lasso = df.coef$gene |> as.character()

set.seed(123)
res.bor = Boruta(as.matrix(df.exp), factor(df.pheno$group), maxRuns=500, doTrace=1)
df.imp = res.bor$ImpHistory |> as.data.frame()
df.plot = reshape2::melt(cbind(x = 1:nrow(df.imp), df.imp), id.var = "x")
df.plot$final = as.character(res.bor$finalDecision)[df.plot$variable]
df.plot$final = ifelse(is.na(df.plot$final), "Shadow", df.plot$final)
df.imp = df.imp[order(apply(df.imp, 2, median), decreasing = F)]
df.plot$variable = factor(df.plot$variable, levels = colnames(df.imp))
df.plot = df.plot[!is.infinite(df.plot$value),]

ggplot(df.plot, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = final), alpha = 0.5) + 
  scale_fill_manual(values = c("Shadow"="blue","Confirmed"="green",
                               "Tentative"="yellow","Rejected"="red"),
                    guide = guide_none()) +
  ylab("Importance") + xlab(NULL) +
  theme_classic2(base_size = 13) + 
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("06.Selection/04.Boruta.Importance.png", width = 6, height = 4, dpi = 300, units = "in", bg = "white")
ggsave("06.Selection/04.Boruta.Importance.pdf", width = 6, height = 4, dpi = 300, units = "in", bg = "white")

gs.bor = res.bor$finalDecision[res.bor$finalDecision == "Confirmed"] |> names()

ggvenn(list(`LASSO` = gs.lasso, `Boruta` = gs.bor),
       stroke_size = 0.5, set_name_size = 5, text_size = 5, show_percentage = F) + 
  ylim(c(-2,2))
ggsave("06.Selection/05.Model.Venn.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("06.Selection/05.Model.Venn.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

gs = intersect(gs.lasso, gs.bor)
write.table(data.frame(gene = gs), "06.Selection/06.Selected.xls", quote = F, sep = "\t", row.names = F)


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,]

df.diff = read.delim2("02.DEG/01.RNAseq.DEG.xls")
df.diff = subset(df.diff, Gene %in% gs)
df.diff$adj.P.Val = as.numeric(df.diff$adj.P.Val)
df.diff$lab = paste0(df.diff$Gene, "\np.adj = ", signif(df.diff$adj.P.Val, 3))
df.diff = df.diff[order(df.diff$adj.P.Val,decreasing = F),]

df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = str_remove_all(df.exp$sample, "[0-9]")
df.exp$group = factor(df.exp$group, levels = c("NC","Disease"))
df.exp$lab = df.diff$lab[match(df.exp$gene, df.diff$Gene)]

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group, color = group), alpha = 0.5) +
  geom_boxplot(aes(color = group), fill = "white", notch = F, width = 0.2) +
  facet_wrap(~lab, scale = "free", ncol = 5) +
  ylab("Expression") + xlab(NULL) +
  ggtitle("RNASeq") +
  scale_fill_jama(guide = guide_none()) + 
  scale_color_jama(guide = guide_none()) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("06.Selection/07.RNAseq.Box.png", width = 7, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("06.Selection/07.RNAseq.Box.pdf", width = 7, height = 4, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,]
df.exp = t(df.exp) |> as.data.frame()

cor.im = cor(df.exp, method = "spearman")
cor.test.im = cor.mtest(df.exp, method = "spearman")$p

df.cor = cbind(gene1 = rownames(cor.im), as.data.frame(cor.im))
df.p = cbind(gene1 = rownames(cor.test.im), as.data.frame(cor.test.im))
df.cor = reshape2::melt(df.cor, id.var = "gene1", variable.name = "gene2", value.name = "cor")
df.p = reshape2::melt(df.p, id.var = "gene1", variable.name = "gene2", value.name = "p")
df.out = merge(df.cor, df.p)
write.table(df.out, "06.Selection/08.RNAseq.Cor.xls", quote = F, sep = "\t", row.names = F)

png("06.Selection/09.RNAseq.Cor.png", width = 6, height = 6, units = "in", bg = "white", res = 300)
par(oma = c(1,1,3,1))
corrplot(cor.im, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = cor.test.im, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
title("RNAseq", outer = T)
dev.off()
pdf("06.Selection/09.RNAseq.Cor.pdf", width = 6, height = 6)
par(oma = c(1,1,3,1))
corrplot(cor.im, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = cor.test.im, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
title("RNAseq", outer = T)
dev.off()


rm(list = ls())
df.exp = read.delim2("cleandata/GSE102673.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
df.pheno = read.delim2("cleandata/GSE102673.pheno.xls")
df.pheno$group = ifelse(df.pheno$tissue.ch1 == "normal skin", "NC", "Disease")
df.exp = df.exp[df.pheno$geo_accession]
gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,]

df.diff = read.delim2("02.DEG/04.GSE102673.DEG.xls")
df.diff = subset(df.diff, Gene %in% gs)
df.diff$P.Value = as.numeric(df.diff$P.Value)
df.diff$lab = paste0(df.diff$Gene, "\np.val = ", signif(df.diff$P.Value, 3))
df.diff = df.diff[order(df.diff$adj.P.Val,decreasing = F),]

df.exp = cbind(gene = rownames(df.exp), df.exp)
df.exp = reshape2::melt(df.exp, id.var = "gene", variable.name = "sample")
df.exp$group = df.pheno$group[match(df.exp$sample, df.pheno$geo_accession)]
df.exp$group = factor(df.exp$group, levels = c("NC","Disease"))
df.exp$lab = df.diff$lab[match(df.exp$gene, df.diff$Gene)]

ggplot(df.exp, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group, color = group), alpha = 0.5) +
  geom_boxplot(aes(color = group), fill = "white", notch = F, width = 0.2) +
  facet_wrap(~lab, scale = "free", ncol = 5) +
  ylab("Expression") + xlab(NULL) +
  ggtitle("GSE102673") +
  scale_fill_jama(guide = guide_none()) + 
  scale_color_jama(guide = guide_none()) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("06.Selection/10.GSE102673.Box.png", width = 7, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("06.Selection/10.GSE102673.Box.pdf", width = 7, height = 4, units = "in", dpi = 300, bg = "white")


rm(list = ls())
df.exp = read.delim2("cleandata/GSE102673.exp.xls", row.names = 1) |> 
  apply(c(1,2), as.numeric) |> as.data.frame()
gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene
df.exp = df.exp[gs,]
df.exp = t(df.exp) |> as.data.frame()

cor.im = cor(df.exp, method = "spearman")
cor.test.im = cor.mtest(df.exp, method = "spearman")$p

df.cor = cbind(gene1 = rownames(cor.im), as.data.frame(cor.im))
df.p = cbind(gene1 = rownames(cor.test.im), as.data.frame(cor.test.im))
df.cor = reshape2::melt(df.cor, id.var = "gene1", variable.name = "gene2", value.name = "cor")
df.p = reshape2::melt(df.p, id.var = "gene1", variable.name = "gene2", value.name = "p")
df.out = merge(df.cor, df.p)
write.table(df.out, "06.Selection/11.GSE102673.Cor.xls", quote = F, sep = "\t", row.names = F)

png("06.Selection/12.GSE102673.Cor.png", width = 6, height = 6, units = "in", bg = "white", res = 300)
par(oma = c(1,1,3,1))
corrplot(cor.im, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = cor.test.im, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
title("GSE102673", outer = T)
dev.off()
pdf("06.Selection/12.GSE102673.Cor.pdf", width = 6, height = 6)
par(oma = c(1,1,3,1))
corrplot(cor.im, method = "circle", diag = T, type = "lower", 
         tl.col = "black", col = colorRampPalette(c("blue","white","red"))(100),
         p.mat = cor.test.im, insig = "blank", addCoef.col = "black",
         number.cex = 0.8)
title("GSE102673", outer = T)
dev.off()


rm(list = ls())
gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c(pal_aaas()(10),"gold4")) +
    #scale_color_aaas() +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank())
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    # scale_color_aaas() +
    scale_color_manual(values = c(pal_aaas()(10),"gold4"))
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcat %in% c("CP:KEGG"))[c(3,4)]

df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> t() |> as.data.frame()
gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene

for(g in gs){
  df.sub = df.exp[-which(colnames(df.exp)==g)]
  exp.g = df.exp[[g]]
  res.cor = cor(df.sub, exp.g) |> as.vector()
  names(res.cor) = colnames(df.sub)
  res.cor = res.cor[order(res.cor, decreasing = T)]
  res.cor = na.omit(res.cor)
  
  res.gsea = GSEA(res.cor, TERM2GENE = df.kegg, pvalueCutoff = 0.05, 
                  seed = 1, pAdjustMethod = "BH", eps = 0)
  top.kegg = res.gsea@result
  write.table(top.kegg, paste0("06.Selection/13.GSEA/", g, ".xls"), quote = F, sep = "\t", row.names = F)
  top.kegg = subset(top.kegg, p.adjust<0.05)
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]
  top.kegg = top.kegg[1:10,]
  top.kegg = na.omit(top.kegg)
  top.kegg = top.kegg[order(top.kegg$NES, decreasing = T),]
  top.kegg = top.kegg$Description
  
  p = gsea.plot(res.gsea, top.kegg, g)
  fn1 = paste0("06.Selection/13.GSEA/", g, ".png")
  fn2 = paste0("06.Selection/13.GSEA/", g, ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
}


rm(list = ls())
df.exp = read.delim2("cleandata/fpkm.xls", row.names = 1) |> 
  apply(c(1,2), \(x) log2(as.numeric(x)+1)) |> as.data.frame()
df.m = msigdbr(category = "H")
gs = split(df.m$gene_symbol, df.m$gs_name)
gsvaPar = ssgseaParam(as.matrix(df.exp), gs)
res.gsva = gsva(gsvaPar)
res.gsva = t(res.gsva) |> as.data.frame()
res.gsva = res.gsva[str_starts(rownames(res.gsva), "Disease"),]

df.exp = t(df.exp) |> as.data.frame()
df.exp = df.exp[rownames(res.gsva),]
df.sel = read.delim2("06.Selection/06.Selected.xls")
df.exp = df.exp[df.sel$gene]

res.cor = linkET::correlate(res.gsva, df.exp, method = "spearman")
df.cor = res.cor$r
df.p = res.cor$p

df.cor.out = cbind(hallmark = rownames(df.cor), as.data.frame(df.cor))
df.p.out = cbind(hallmark = rownames(df.p), as.data.frame(df.p))
df.cor.out = reshape2::melt(df.cor.out, id.var = "hallmark", variable.name = "gene", value.name = "cor")
df.p.out = reshape2::melt(df.p.out, id.var = "hallmark", variable.name = "gene", value.name = "p")
df.out = merge(df.cor.out, df.p.out)
write.table(df.out, "06.Selection/14.HALLMARK.Cor.xls", quote = F, sep = "\t", row.names = F)

df.out = subset(df.out, p < 0.05)
ggplot(df.out, aes(x = gene, y = hallmark)) + 
  geom_point(aes(color = cor, size = abs(cor))) + 
  geom_text(aes(label = signif(cor, 2)), size = 3) +
  scale_color_gradient2(low = "blue", high = "red", limits = c(-1,1)) + 
  xlab(NULL) + ylab(NULL) +
  scale_size_continuous(range = c(0,9), guide = guide_none(), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1))
ggsave("06.Selection/15.HALLMARK.Cor.png", width = 7, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("06.Selection/15.HALLMARK.Cor.pdf", width = 7, height = 6, units = "in", dpi = 300, bg = "white")




