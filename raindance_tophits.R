load(file.path(Sys.getenv("HOME"), "Downloads/raindance/2015-07-23_processed.Rdata"))

setwd(file.path(Sys.getenv("HOME"), "git/raindance"))
source("raindance_tophits_functions.R")


vcf$key <- get_key(vcf)
vcfout$key <- get_key(vcfout)
vcf$sitekey <- get_site_key(vcf)

vcf$P50[which(vcf$P50 == "1")] <- NA

vcf.rep <- get_replicated_sites(vcf)

vcf$Rep <- "Uni"
vcf$Rep[which(vcf$sitekey %in% vcf.rep$sitekey)] <- "Rep"

table(vcf$Rep)

repliatesonly <- vcf[which(vcf$Rep == "Rep"),]

vcf <- cbind(indc, libc, vcf)

vcf <- move_sample_cols_to_last(vcf)
vcf <- convert_delimiter_to_semicolon(vcf)

cat(dim(vcf), "\n")

write.table(
    vcf,
    file=paste0("postqc.", nrow(vcf), "list.", timestamp, ".tsv"),
    sep="\t", col.names=TRUE, row.names=FALSE)


#########################################################
## count the number of calls
#########################################################
ncallp20 <- sum(sapply(vcf$P20, function(x){
    length(unlist(strsplit(x, ";")))
}))
cat("number of calls at P20:", ncallp20, "\n")


get_num_of_calls <- function(invcf, pcol){
    ncallpcall <- 0
    for(i in 1:nrow(invcf)){
        pcall <- invcf[i, pcol]
        if(is.na(pcall)) {
            next
        } else {
            ncallpcall <- ncallpcall + length(unlist(strsplit(pcall, ";")))
        }
    }
    cat("number of calls at ", pcol,":", ncallpcall, "\n")
    ncallpcall
}

ncallp20 <- get_num_of_calls(vcf, "P20")
ncallp50 <- get_num_of_calls(vcf, "P50")

#########################################################
## top list in response to the following

## 1. Starting with your post-blacklisting list (6054 calls), apply
## the following filters

## Remove variants 1bp after primer
## VAF upper limit of 3%
## Remove site if median coverage <5000

## Send this new list to us


## Apply final filtering criteria to get to 413 list

## Select variants called in more than one sample and/or more than one
## amplicon

## Exclude variant if all calls are made from 1 library and the number
## of calls is 3 or more (unless it is called in >1 amplicon)

## I think it’s important to get an idea how much the stage 1
## filtering excludes variants and how much the stage 2 does.

## As per original email, please can you alter the P20 and P50 columns
## to avoid the comma issues, and if possible to have columns with
## biopsy name.
#########################################################

## timestamp <- format(Sys.time(), "%Y_%b_%d_%H-%M-%S")
timestamp <- format(Sys.time(), "%Y_%b_%d_%H")

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

vcf.stage1 <- move_sample_cols_to_last(vcf[iselect.stage1, ])
vcf.stage1 <- convert_delimiter_to_semicolon(vcf.stage1)

cat(dim(vcf.stage1), "\n")

write.table(
    vcf.stage1,
    file=paste0("topsites.stage1.", nrow(vcf.stage1), "list.", timestamp, ".tsv"),
    sep="\t", col.names=TRUE, row.names=FALSE)

########################
## stage2

## Select variants called in more than one sample and/or more than one
## amplicon

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

vcf.stage2 <- move_sample_cols_to_last(vcf[iselect.stage2, ])
vcf.stage2 <- convert_delimiter_to_semicolon(vcf.stage2)

## just check if the sites of interest are in the list
stopifnot(
    nrow(vcf.stage2[which(vcf.stage2$Pos == 112888199 & vcf.stage2$Alt == "A"),
                   metric_cols(vcf.stage2)]) == 2)

stopifnot(check_no_singletons(vcf.stage2))

stopifnot(check_know_replicates(vcf.stage2))

write.table(
    vcf.stage2,
    file=paste0("topsites.stage2.", nrow(vcf.stage2), "list.", timestamp, ".tsv"),
    sep="\t", col.names=TRUE, row.names=FALSE)

#########################################################
#########################################################

## save.image(paste0("manuscript.", timestamp, ".Rdata"))

save.image("manuscript.2018_Mar_03_21.Rdata")

stop("done")





query.ret <- vcf[which(vcf$Amp == "64:RET_5:chr10:43614940-43615035"),]
query.ret <- query.ret[which(query.ret$Pos == 43615026),]

query.nras <- vcf[which(vcf$Amp == "3:NRAS_3:chr1:115251240-115251336"),]
query.nras <- query.ret[which(query.ret$Pos == 115251320),]


## why they didn't make it into the final 400ish list

query <- read.csv("RD_6054_exclusioncriteria_medcoverage.csv", stringsAsFactors=FALSE)

query <- query[, c( "X.2.POS", "X.6.Amp.1", "X.4.REF", "X.5.ALT")]
colnames(query) <- c("Pos", "Amp", "Ref", "Alt")

query.dat <- merge(query,
		   vcf[, c("Pos", "Amp", "Ref", "Alt",
			   "Rep", "PRank", "MdCov", "RelPos",
			   "P20LibMax", "P20SizeLib",
			   "P20SizeInd", "MxRate", "MxRateRank")],
		   by.x=c("Pos", "Amp", "Ref", "Alt"),
		   by.y=c("Pos", "Amp", "Ref", "Alt"), all.x=TRUE, sort=FALSE)

query.dat$"BadLibrary" <- "No"
query.dat.badlibidx <- which(query.dat$P20LibMax>=3 & query.dat$P20SizeLib == 1)
query.dat$"BadLibrary"[query.dat.badlibidx] <- "Yes"

idx <- which(query.dat$MdCov>5000 & query.dat$MxRate>0 & query.dat$MxRate<0.03 & query.dat$RelPos>2 & query.dat$MxRateRank < 1800 & (query.dat$Rep == "Rep" | (query.dat$PRank < 1500)))
idx <- idx[!idx %in% query.dat.badlibidx]

query.dat$"Included_in_463list" <- "No"
query.dat$"Included_in_463list"[idx] <- "Yes"

write.table(query.dat, "query.36sites.dat.tsv", col.names=TRUE, row.names=FALSE, sep="\t")

toplist <- read.csv("anno.rank.top.txt", sep="\t")


query.dat <- merge(query.dat,
		   toplist[, c("Pos", "Amp", "Ref", "Alt")],
		   by.x=c("Pos", "Amp", "Ref", "Alt"),
		   by.y=c("Pos", "Amp", "Ref", "Alt"), all.x=TRUE, sort=FALSE)



####################################################
## Make the top list
####################################################

badlibidx <- which(vcf$P20LibMax>=3 & vcf$P20SizeLib == 1)

## not clear what criteria were used
topidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & (vcf$Rep == "Rep" | (vcf$PRank < 1500 & (vcf$P20SizeInd >=2))))
topidx <- topidx[!topidx %in% badlibidx]
length(topidx)

## This is the criteria learnt from Hannah's list
topidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$MxRate<0.03 & vcf$RelPos>2 & vcf$MxRateRank < 1800 & (vcf$Rep == "Rep" | (vcf$PRank < 1500)))
topidx <- topidx[!topidx %in% badlibidx]
length(topidx)

vcf$top <- "No"
vcf$top[topidx] <- "Yes"

vcfout$top <- "No"
vcfout$top[topidx] <- "Yes"

write.table(vcfout[topidx,], file="463list.tsv", sep="\t", col.names=TRUE, row.names=FALSE)

## Rebuild toplist
toplist.redo <- vcf[topidx,]
toplist.redo$key <- get_key(toplist.redo)

cat("size of the rebuilt toplist ", dim(toplist.redo), "\n")

## Hannah's toplist
toplist.hr <- read.csv("RDanno.priority_24March15_Oct15_edit.csv", stringsAsFactors=FALSE)
toplist.hr$key <- get_key(toplist.hr)
cat("size of Hannah's toplist ", dim(toplist.hr), "\n")

## Hannah's toplist
toplist.anne <- read.csv("RDanno.rank.top_10March15.csv", stringsAsFactors=FALSE)
toplist.anne$key <- get_key(toplist.anne)
cat("size of Anne's toplist ", dim(toplist.anne), "\n")

####################################################
## Reverse engineer from Hannah's list
####################################################

table(vcf[which(vcf$key %in% toplist.anne$key), "top"])

summary(vcf[which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$MxRate<0.03 & vcf$RelPos>2 & vcf$Rep != "Rep"), metric_cols(vcf)])

summary(toplist.anne[, metric_cols(toplist.anne)])

toplist.hr[, metric_cols(toplist.hr)]

####################################################
## Compare the two toplists
## one from Hannah and one regenerated in Feb 2018
## both have 463 entries
####################################################

notfound <- setdiff(toplist.hr$key, toplist.redo$key)

notfound[grep("FGFR2", notfound)]

vcf[which(vcf$Pos == "120458776"), metric_cols(vcf)]

setdiff(toplist.redo$key, toplist.hr$key)

####################################################
## Top list in the team drive at WIMM
## /Volumes/Wilkie/WilkieGr/Sperm/Raindance/anno.rank.top.txt
## renamed to /Users/zhihao/Downloads/raindance/anno.rank.top.txt_wimm
####################################################

toplist.wimm <- read.csv("/Users/zhihao/Downloads/raindance/anno.rank.top.txt_wimm", sep="\t")

toplist.wimm$key <- get_key(toplist.wimm)
[
which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & (vcf$Rep == "Rep" | (vcf$PRank < 1500 & (vcf$P20SizeInd >=2))))

summary(toplist.wimm[, metric_cols(toplist.wimm)])

setdiff(toplist.wimm$key, toplist.hr$key)


####################################################
## Notch3 example
## the current list has one more site with a different
## mutation than the previous list
####################################################

vcf[which(vcf$Pos == 15303217), c(metric_cols(vcf), which(colnames(vcf) %in% c(52,55,56,69,62,65)))]

toplist.hr[which(toplist.hr$Pos == 15303217), metric_cols(toplist.hr)]
toplist.redo[which(toplist.redo$Pos == 15303217), metric_cols(toplist.redo)]

write.table(toplist.redo[which(toplist.redo$Pos == 15303217), c(metric_cols(toplist.redo), 99)],
	    file="Nortch3.example.tsv", sep="\t", col.names=TRUE, row.names=FALSE)

amp <- read.csv("/Users/zhihao/Downloads/raindance_data/489:Notch3_5:chr19:15303178-15303273/489:Notch3_5:chr19:15303178-15303273.ncc.csv")
amp[which(amp$Pos==15303217), 20]

amp.t <- read.csv("/Users/zhihao/Downloads/raindance_data/489:Notch3_5:chr19:15303178-15303273/489:Notch3_5:chr19:15303178-15303273.T.csv")
amp.t[which(amp.t$Pos==15303217),]

amp.dat <- read.csv("/Users/zhihao/Downloads/raindance_data/489:Notch3_5:chr19:15303178-15303273/489:Notch3_5:chr19:15303178-15303273.calls.qtsd.csv", header=TRUE)

amp.dat[which(amp.dat$Pos == 15303217),]

####################################################
## LRP5 example
## 141:LRP5_17:chr11:68153947-68154100
## the current list has one more site with a different
## mutation than the previous list
####################################################

vcf[which(vcf$Pos == 68154097), c(metric_cols(vcf), which(colnames(vcf) %in% c(173,175,181,200,206,209,224)))]

load("/Users/zhihao/Downloads/raindance/example_data/141:LRP5_17:chr11:68153947-68154100/141:LRP5_17:chr11:68153947-68154100.Rdata")

write.table(vcf[which(vcf$Pos == 68154097), metric_cols(vcf)],
	    file="LRP5.example.tsv", sep="\t", col.names=TRUE, row.names=FALSE)
