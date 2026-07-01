rm(list = ls())
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggradar)

gs = read.delim2("06.Selection/06.Selected.xls")
gs = gs$gene

for(g in gs){
  fn = paste0("interm/ctd/", g, ".csv")
  df = read.csv(fn)
  df = df[c(1,5)]
  df = df[!duplicated(df),]
  df = df[order(df$Inference.Score, decreasing = T),]
  df = df[1:15,]
  df$Disease.Name = factor(df$Disease.Name, levels = rev(df$Disease.Name))
  lim.min = round(round(min(df$Inference.Score),0)-5,-1)
  lim.max = round(round(max(df$Inference.Score),0)+15,-1)-5
  if(lim.max > 300) lim.max = lim.max + 50
  
  ggplot(df, aes(x = Inference.Score, y = Disease.Name)) + 
    geom_segment(aes(yend = Disease.Name, xend = lim.min)) +
    geom_text(aes(label = round(Inference.Score,2)), hjust = -0.3) +
    geom_point(aes(size = Inference.Score, color = Inference.Score)) +
    xlab("Inference Score") + 
    ylab("Disease") + ggtitle(g) +
    scale_x_continuous(expand = c(0,0), limits = c(lim.min, lim.max)) +
    scale_size_continuous(range = c(3,6), guide = guide_none()) +
    scale_color_gradient2(low = "yellow",high = "darkred", guide = guide_none(), limits = c(lim.min, lim.max)) +
    theme_classic2() + 
    theme(plot.title = element_text(hjust = 0.5))
  fn1 = paste0("10.Disease/", g, ".png")
  fn2 = paste0("10.Disease/", g, ".pdf")
  ggsave(fn1, width = 8, height = 6, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, width = 8, height = 6, units = "in", dpi = 300, bg = "white")
}
