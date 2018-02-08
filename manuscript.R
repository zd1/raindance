library(data.table)
library(compare)

extlist <- data.frame(fread("/home/zhihao/Downloads/All_RDcall_18Feb15.csv"))

curated <- data.frame(fread("/home/zhihao/Downloads/anno.rank.txt", header=TRUE))

comparison <- compare(extlist[, c("Chr", "Start", "End", "Ref", "Alt")],
                      curated[, c("Chr", "Start", "End", "Ref", "Alt")],
                      allowAll=TRUE)
