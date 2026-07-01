rm(list = ls())

df.mirna = read.csv("09.Regulation/01.miRNA.csv")
cnt = table(df.mirna$ID)
cnt = cnt[cnt == 5]
df.mirna = df.mirna[df.mirna$ID %in% names(cnt),]

df.lnc = read.csv("09.Regulation/02.lncRNA.csv")
df.lnc = df.lnc[!duplicated(df.lnc),]
df.lnc = df.lnc[df.lnc$ID %in% df.mirna$ID,]
cnt = table(df.lnc$Target)
cnt = cnt[cnt == 5]
df.lnc = df.lnc[df.lnc$Target %in% names(cnt),]

df.mirna = df.mirna[c(1,3)]
df.lnc = df.lnc[c(1,3)]

df.full = rbind(df.mirna, df.lnc)
write.table(df.full, "interm/mirna.net.txt", quote = F, sep = "\t", row.names = F)
