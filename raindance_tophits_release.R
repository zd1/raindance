vcf <- read.csv("postqc.6054list.2018_Mar_03_21.tsv", sep="\t", header=TRUE)

########################
## stage1
## Remove variants 1bp after primer
## VAF upper limit of 3%
## Remove site if median coverage <5000
########################

iselect.stage1 <- which(
    vcf$MdCov>5000 &
    vcf$RelPos>1 &
    vcf$MxRate<0.03)
iselect.stage1 <- setdiff(iselect.stage1, badlibidx)
cat("length of selection in stage 1:", length(iselect.stage1), "\n")

vcf.stage1 <- vcf[iselect.stage1, ]

cat(dim(vcf.stage1), "\n")

########################
## stage2

## Select variants called in more than one sample and/or more than one amplicon
## Exclude variant if all calls are made from 1 library and the number
## of calls is 3 or more (unless it is called in >1 amplicon)
########################

badlibidx <- which(vcf$P20LibMax>=3 & vcf$P20SizeLib == 1)
badlibidx <- setdiff(badlibidx, which(vcf$Rep == "Rep"))
iselect.stage2 <- which(
    vcf$MdCov>5000 &
    vcf$RelPos>1 &
    vcf$MxRate<0.03 & (vcf$Rep == "Rep" | vcf$P20SizeInd >=2))
iselect.stage2 <- setdiff(iselect.stage2, badlibidx)
cat("length of selection:", length(iselect.stage2), "\n")

vcf.stage2 <- vcf[iselect.stage2, ]

save.image("manuscript.2018_Mar_03_21.Rdata")
