rm(list = ls())
library(tidyverse)
library(ggsci)
library(ggpubr)

df.loc = read.csv("interm/tableExport.csv")
df.loc = df.loc[-c(1,3)]
df.loc = reshape2::melt(df.loc, id.var = "FastaID")

ggplot(df.loc, aes(x = FastaID, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_nejm(name = "Location") +
  scale_y_continuous(expand = c(0,0)) +
  xlab(NULL) + ylab("Proportion") +
  theme_classic2()
ggsave("08.Location/01.Cell.Location.png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("08.Location/01.Cell.Location.pdf", width = 7, height = 5, units = "in", dpi = 300, bg = "white")


rm(list = ls())
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=3, tracks.outside=0)

df.gtf = read.delim2("../database/Human.GRCh38.p13.annot.tsv.gz")
hub = read.delim2("06.Selection/06.Selected.xls")
df.gtf = df.gtf[df.gtf$Symbol %in% hub$gene,]
df.gtf$Chromosome = str_remove(df.gtf$ChrAcc, "NC_0*") |> str_remove("\\..*")
df.gtf$Chromosome = paste0("chr", df.gtf$Chromosome)
df.gtf = df.gtf[c("Chromosome","ChrStart","ChrStop","Symbol")]
colnames(df.gtf) = c("Chromosome","chromStart","chromEnd","Gene")
df.gtf$chromStart = as.numeric(df.gtf$chromStart)
df.gtf$chromEnd = as.numeric(df.gtf$chromEnd)

png("08.Location/02.Genome.Location.png", width = 5, height = 5, units = "in", res = 300, bg = "white")
params = RCircos.Get.Plot.Parameters()
params$text.size = 0.5
RCircos.Reset.Plot.Parameters(params)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(df.gtf, 1, "in")
RCircos.Gene.Name.Plot(df.gtf, 4, 2, "in")
dev.off()

pdf("08.Location/02.Genome.Location.pdf", width = 5, height = 5)
params = RCircos.Get.Plot.Parameters()
params$text.size = 0.5
RCircos.Reset.Plot.Parameters(params)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(df.gtf, 1, "in")
RCircos.Gene.Name.Plot(df.gtf, 4, 2, "in")
dev.off()

