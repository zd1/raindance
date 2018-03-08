
library(ggplot2)
library(reshape)

wkdir="/Users/zd1/cloud/data/raindance/pileup/sum/calls/annovar/release_2015-07-23/"
setwd(wkdir)
samples <- read.table("/Users/zd1/cloud/data/raindance/meta/samples.csv", sep=",", header=TRUE)

timestamp <- format(Sys.time(), "%Y_%b_%d_%H-%M-%S")
save.image(paste0(timestamp, ".Rdata"))
######################################################################################
## Unannotated VCF table made by bcftools, we filtered on Q20, so have 6k variants
######################################################################################

# vcf
vcf=read.table("/Users/zd1/volumn/wt/raindance/call/simon/all.flt.anvinput", header=TRUE, sep="\t", comment.char="", stringsAsFactors=FALSE)
colnames(vcf)[1:10] <- c("Chrm", "Pos", "Pos", "Ref", "Alt",
                         "Amp", "MaxP", "P20", "P50", "X")
colnames(vcf)[11:298] <- seq(1,288)

## Annotation for the same set of variants
anno=read.table("/Users/zd1/volumn/wt/raindance/call/simon/all.flt.anvinput.annovar.hg19_multianno.txt", header=TRUE, sep="\t")
anno.vcf <- cbind(anno, vcf)

## Amplicon
ampcat=read.table("/Users/zd1/cloud/data/raindance/meta/ampcat.csv", header=FALSE, sep=",")
colnames(ampcat) <- c("AmpID", "Amp", "CatNum")
amp <- read.table("/Users/zd1/cloud/data/raindance/meta/amplicons.csv", header=TRUE, sep=",")
amp$ts <- sapply(1:dim(amp)[1], function(i){
    as.numeric(paste(strsplit(as.character(amp$Target.Start)[i],",")[[1]],collapse=""))
           })
amp$te <- sapply(1:dim(amp)[1], function(i){
    as.numeric(paste(strsplit(as.character(amp$Target.End)[i],",")[[1]],collapse=""))
           })

amp$Amp <- sapply(1:dim(amp)[1], function(i){
    s <- as.numeric(paste(strsplit(as.character(amp$Chr.Start)[1],",")[[1]],collapse=""))
    e <- as.numeric(paste(strsplit(as.character(amp$Chr.End)[1],",")[[1]],collapse=""))
    chrm <- paste("chr",as.character(amp$Chr[1]), sep="")
    aid <- as.character(amp$Amplicon.ID[i])
    ex <- as.character(amp$Customer.TargetID[i])
    paste(aid,":", ex, ":",chrm,":",s,"-",e,sep="")
})

save.image("2015-07-23.Rdata")
############################################################################
############################################################################
get_key <- function(dat){
    apply(dat[, c("Chrm", "Pos", "Amp", "Ref", "Alt")], 1, function(x){
	paste0(trimws(x), collapse="_")
    })
    ## apply(dat[, c("Chrm", "Pos", "Amp", "Ref")], 1, function(x){
    ##     paste0(trimws(x), collapse="_")
    ## })
}

get_site_key <- function(dat){
    apply(dat[, c("Chrm", "Pos", "Ref", "Alt")], 1, function(x){
	paste0(trimws(x), collapse="_")
    })
    ## apply(dat[, c("Chrm", "Pos", "Amp", "Ref")], 1, function(x){
    ##     paste0(trimws(x), collapse="_")
    ## })
}

convert_delimiter_to_semicolon <- function(invcf){
    invcf$P20 <- sapply(invcf$P20, function(x){gsub(",", ";", x)})
    invcf$P50 <- sapply(invcf$P50, function(x){gsub(",", ";", x)})
    invcf
}

############################################################################
## Mutation Pattern
############################################################################
## load("2015-07-23.Rdata")
vcf <- anno.vcf

a.let<-c("A", "C", "G", "T");
mut <- array(c(1,2,1,3,1,4,2,1,2,3,2,4,3,1,3,2,3,4,4,1,4,2,4,3), c(2,12))
mutptn <- c()
for(i in 1:12){
    from = mut[1, i]
    to = mut[2, i]
    from2=mut[1, 13-i]
    to2=mut[2, 13-i]
    if(i<=6){
        mutptn <- rbind(mutptn, c(from, to, a.let[from], a.let[to],
                              paste(a.let[from],">", a.let[to],
                                    "(", a.let[from2], ">", a.let[to2], ")", sep="")))
    }else{
        mutptn <- rbind(mutptn, c(from, to, a.let[from], a.let[to],
                              paste(a.let[from2],">", a.let[to2],
                                    "(", a.let[from], ">", a.let[to], ")", sep="")))
    }
}

varmut <- apply(vcf, 1, function(x){
     mutptn[which(as.character((x[4])) ==  mutptn[,3] & as.character((x[5])) ==  mutptn[,4]),5]
})


vcf$Mut <- varmut


############################################################################
## Expected number of variants by chance according to amp length
############################################################################

amplen <- sapply(1:dim(vcf)[1], function(x){
    y <- as.integer(strsplit(strsplit(as.character(vcf$Amp)[x], ":")[[1]][4],"-")[[1]])
    abs(y[1] - y[2])
})
vcf$Amplen <- amplen

amplen <- sapply(1:dim(vcf)[1], function(x){
    y <- as.integer(strsplit(strsplit(as.character(vcf$Amp)[x], ":")[[1]][4],"-")[[1]])
    abs(y[1] - y[2])
})
vcf$Amplen <- amplen


############################################################
## Add categories
############################################################

catname <- c("1","FiveGenes", "2", "MAPK", "3", "COSMIC", "4", "Disease", "c", "Control")
catname <- data.frame(t(array(catname,c(2,5))))
colnames(catname) <- c("CatNum", "Cat")

ampcat$Amplen  <- sapply(1:dim(ampcat)[1], function(x){
    y <- as.integer(strsplit(strsplit(as.character(ampcat$Amp)[x], ":")[[1]][4],"-")[[1]])
    abs(y[1] - y[2])
})

ampcat$AmplenRatio <- ampcat$Amplen/sum(ampcat$Amplen)

ampcat$Cat <- catname$Cat[match(as.character(ampcat$CatNum), as.character(catname$CatNum))]
vcf$Cat <- ampcat$Cat[match(as.character(vcf$Amp), as.character(ampcat$Amp))]
vcf$AmplenRatio <- ampcat$AmplenRatio[match(as.character(vcf$Amp), as.character(ampcat$Amp))]

vcf$Cat <- ampcat$Cat[match(as.character(vcf$Amp), as.character(ampcat$Amp))]
vcf$AmplenRatio <- ampcat$AmplenRatio[match(as.character(vcf$Amp), as.character(ampcat$Amp))]

pdf("freq_p20.pdf")
barplot(table(vcf$Cat), main="Freq at P20")
dev.off()

relpos <- sapply(1:dim(vcf)[1], function(i){
    min(abs(as.integer(vcf$Pos[i]) - as.integer(strsplit(strsplit(as.character(vcf$Amp)[i], ":")[[1]][4],"-")[[1]])))
})

vcf$RelPos <- relpos
#vcfout$RelPos <- relpos

################################################################
## Add Median Coveage, Max allele frequency, whether replicated
################################################################

altcol <- 5
sampstart <- 15
n.sam <- 288

sami <- which(colnames(vcf) %in% seq(1,288))

vcfmeta <- apply(vcf, 1, function(x){
    x <- x[sami] # samples
    vsc <- lapply(x, function(y) {as.numeric(strsplit(strsplit(as.character(y[1]), ":")[[1]][1], ",")[[1]])})
    mdcov <- median(sapply(vsc, sum))
    mxrate <- max(sapply(vsc, function(x){xi <- order(x, decreasing=TRUE)[2]; x[xi]/sum(x)}), na.rm=TRUE)
    c(ref, alt, mdcov, mxrate)
})

## vcf$Ref1 <- vcfmeta[1,]
## vcf$Alt1 <- vcfmeta[2,]
vcf$MdCov <- as.numeric(vcfmeta[3,])
vcf$MxRate <- as.numeric(vcfmeta[4,])
vcf$MxRateRank <- rank(-vcf$MxRate)
vcf$MaxP[which(!is.finite(vcf$MaxP))] <- 300
vcf$PRank <- rank(-vcf$MaxP)

vcf$PvlRank <- ave(-log10(vcf$MaxP), vcf$Cat, FUN=rank)
vcf$logRate <- log10(vcf$MxRate)
vcf$RateRank <- ave(-vcf$logRate, vcf$Cat, FUN=rank)

vcf$key <- get_key(vcf)
vcf$sitekey <- get_site_key(vcf)

# Replicated - a variant identified in overlapping amplicons in *same* samples
vcf.rep <- vcf[,c("Chrm", "Pos", "Alt", "sitekey", "P20")]
s <- strsplit(as.character(vcf.rep$P20), split = ",")
groupsizes <- sapply(s, length)
vcf.rep.exd <- data.frame(Chrm=rep(vcf.rep$Chrm, groupsizes),
                          Pos=rep(vcf.rep$Pos, groupsizes),
                          Alt=rep(vcf.rep$Alt, groupsizes),
                          sitekey=rep(vcf.rep$sitekey, groupsizes),
                          P20 = unlist(s))

repi.exd <- which(duplicated(vcf.rep.exd) | duplicated(vcf.rep.exd, fromLast = TRUE))
vcf.rep <- vcf.rep.exd[repi.exd,]
vcf.rep <- unique(vcf.rep[, c("Chrm", "Pos", "Alt", "sitekey")])

known.singletons <- data.frame(do.call("rbind", strsplit(c("chr11_113934509_C", "chr11_113934514_C", "chr11_113934525_C", "chr11_113934547_A", "chr11_113934548_G", "chr11_113934563_C", "chr11_113934581_C", "chr11_113934606_G", "chr11_113934631_T", "chr11_113934640_C", "chr11_113934667_C", "chr11_113934689_C", "chr11_113934718_A", "chr11_113934734_C", "chr11_113934739_G", "chr11_113934697_T", "chr11_113934706_C", "chr11_113934728_G", "chr11_113934732_C", "chr11_113934738_A", "chr11_113934752_C", "chr11_113934758_C", "chr11_113934760_T", "chr11_113934762_C", "chr11_113934776_T", "chr11_113934777_A", "chr11_113934782_A", "chr11_113934786_G", "chr11_113934801_G", "chr11_113934806_A", "chr11_113934807_A"), split="_")))
colnames(known.singletons) <- c("Chrm", "Pos", "Alt")

stopifnot(
    nrow(merge(known.singletons,
               vcf.rep,
               by.x=c("Chrm", "Pos", "Alt"),
               by.y=c("Chrm", "Pos", "Alt"))) == 0)

vcf$Rep <- "Uni"
vcf$Rep[which(vcf$sitekey %in% vcf.rep$sitekey)] <- "Rep"

stopifnot(
    nrow(merge(known.singletons,
               vcf[which(vcf$Rep == "Rep"), c("Chrm", "Pos", "Alt")],
               by.x=c("Chrm", "Pos", "Alt"),
               by.y=c("Chrm", "Pos", "Alt"))) == 0)
table(vcf$Rep)

repliatesonly <- vcf[which(vcf$Rep == "Rep"),]

write.table(
    repliatesonly,
    file=paste0("replicates.", nrow(repliatesonly), "list.", timestamp, ".tsv"),
    sep="\t", col.names=TRUE, row.names=FALSE)



save.image("2015-07-23_processed.Rdata")


# libraries
libs <- unique(samples$Illumina.lane)
libc <- array(0, c(dim(vcf)[1], length(libs)))

# individuals
inds <- unique(samples$individual)
indc <- array(0,c(dim(vcf)[1],length(inds)))

# add counts
for(i in seq(1,dim(vcf)[1])){
    sam_i <- samples[match(as.integer(strsplit(as.character(vcf$P20[i]), ",")[[1]]), samples$sid),]
    inds_i <- sam_i$individual
    lib_i <- sam_i$Illumina.lane
    for(j in match(inds_i,inds)){
        indc[i,j] = indc[i,j] + 1
    }
    for(k in match(lib_i, libs)){
        libc[i,k] = libc[i,k] + 1
    }
}

colnames(indc) <- paste("Ind",inds,sep="")
colnames(libc) <- libs

indc <- data.frame(indc)
libc <- data.frame(libc)

vcf$P20SizeInd <- apply(indc,1,max)
vcf$P20SizeLib <- apply(libc, 1, function(x){length(x[which(x>0)])})
vcf$P20LibMax <- apply(libc,1,max)

## add mut patterns
varmut <- apply(vcf[,c("Ref", "Alt")], 1, function(x){
     mutptn[which(as.character((x[1])) ==  mutptn[,3] & as.character((x[2])) ==  mutptn[,4]),5]
})

vcf$Mut <- varmut

## output table
features <- colnames(vcf)[which(! colnames(vcf) %in% seq(1,288))]
features <- features[!features %in% c("Chr", "Start", "End", "X",
                                       "Chrm", "Pos", "Ref1", "Alt1")]
features <- c(features, seq(1,288))

vcfout <- cbind(vcf[, c("Chrm", "Pos", "Ref", "Alt")],
                indc, libc, vcf[,features])

write.table(vcfout, "vcfout.csv", col.names=TRUE, row.names=FALSE, sep=",")

save.image("2015-07-23_processed.Rdata")

############################################################
## Mut Pattern Per lane
############################################################

vcfmutptn <- vcfout[, grep("WTCHG",colnames(vcfout))]
vcfmutptn$Mut <- vcfout$Mut

vcfmutptnlib <- melt(vcfmutptn, id=c("Mut"))
colnames(vcfmutptnlib)[2:3] <- c("Lib", "Count")
vcfmutptnlib$Count <- as.numeric(vcfmutptnlib$Count)
vcfmutptnlib$Mut <- as.factor(vcfmutptnlib$Mut)

agg <- aggregate(vcfmutptnlib$Count, by=list(vcfmutptnlib$Lib, vcfmutptnlib$Mut), FUN=sum)
colnames(agg) <- c("Lib", "Mut", "Count")

pdf("lib.mutptn.pdf")
ggplot(agg, aes(x=Lib, y=Count, fill=Mut)) + geom_bar(stat="identity") +
            ggtitle("Type of Mutations per Library") +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1)),
                      axis.title.x = element_text(size = rel(1.8)),
                      axis.title.y = element_text(size = rel(1.8)))
dev.off()


vcfmutptnlibfrq <- data.frame(t(apply(vcfmutptn[,1:16], 1, function(x){x/sum(x)})))
vcfmutptnlibfrq$Mut <- vcfmutptn$Mut
vcfmutptnlibfrq <- melt(vcfmutptnlibfrq, id=c("Mut"))
vcfmutptnlib$frq <- vcfmutptnlibfrq$value

vcfmutptnlibmelt <- melt(vcfmutptnlib, id=c("Mut", "Lib"))


pdf("lib.mutptn.samplecount.pdf")
ggplot(vcfmutptnlibmelt, aes(x=Lib, y=value, group=Mut, color=Mut)) + geom_point() +
    facet_grid(variable ~ .,  scales = "free") +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1)),
                      axis.title.x = element_text(size = rel(1.8)),
                      axis.title.y = element_text(size = rel(1.8)))
dev.off()

############################################################
## Pick Sites
############################################################
load("2015-07-23_processed.Rdata")

## fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & (vcf$PvlRank<100 | vcf$RateRank<100 | vcf$Rep == "Rep") &  vcf$P20Size>1)


badlibidx <- which(vcf$P20LibMax>=3 & vcf$P20SizeLib == 1)

topidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & vcf$MxRate <0.03 &  vcf$MxRateRank < 3000 & (vcf$Rep == "Rep" | vcf$PRank< 5000))
topidx <- topidx[!topidx %in% badlibidx]
vcftop <- vcf[topidx,]
length(unique(vcftop[, "Pos"]))

vcftop[which(vcftop$Pos==15303236), c("MaxP", "MxRate", "MxRateRank")]

write.table(vcfout[topidx, ], "vcfout.csv", col.names=TRUE, row.names=FALSE, sep=",")

checkvariants <- c(123276893, 30729984, 1807371, 140453132, 113934759, 39249807)
checkvariants %in% vcf$Pos[topidx]



#####################################################################
## Make a file for Simon in a format similar to the CGP eyelid paper
####################################################################

res <- c()
p20i <- which(colnames(vcf) == "P20")

for(i in topidx){
    x <- vcf[i,]
    sami <- strsplit(lapply(x[p20i], as.character)[[1]], ",")[[1]]
    for(s in sami){
        d <- strsplit(lapply(x[which(colnames(vcf) == s)], as.character)[[1]], ":")[[1]]
        ct <- as.numeric(strsplit(d[1], ",")[[1]])
        pv <- as.numeric(strsplit(d[2], ",")[[1]])
        cv <- sum(ct)
        mt <- ct[which(x$Alt == a.let)]
        p <- pv[which(x$Alt == a.let)]
        row <- c(s, x$Chrm, x$Pos, as.character(x$Ref),
                 as.character(x$Alt), mt, cv, mt/cv,
                 as.character(p),
                 as.character(x$MaxP),
                 as.character(x$Amp),
                 as.character(x$ExonicFunc.ensGene),
                 as.character(x$AAChange.refGene))
        res <- rbind(res, row)
    }
}

resd <- as.data.frame(res)

fields <- c("sid", "Chrm", "Pos", "Ref", "Alt",
                   "MutCount", "Coverage", "Rate", "mlogP", "MaxP", "Amplicon",
                    "ExonicFunc", "AminoAcidChange")
colnames(resd) <- fields

resd.sam <- merge(resd, samples, by.x=c("sid"), by.y=c("sid"))
resd.out <- resd.sam[,c(fields, "genomes.per.bit", "Ind.Slide.sample")]

write.table(resd.out, "variants.csv", col.names=TRUE, row.names=FALSE, sep=",")

save.image("2015-07-23_processed.Rdata")



############################################################
## Enrichment
############################################################
expcat <- countenrich(vcf[fltidx,])
pdf(paste("expcounts_top.",n,".pdf",sep=""))
ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("MinCov>5000 & (PvlRank<100 or RateRank<100 or Rep) & Size>1"," n=",n, sep=""))
dev.off()

countenrich <- function(vcf){
    n.var <- dim(vcf)[1]
    cats <- levels(vcf[,"Cat"])
    expcat <- c()
    #totallen <-  sum(unique(vcf[, c("Amp","Amplen")])[2])
    for(c in cats){
        ## catlen <- (unique(vcf[which(as.character(vcf$Cat) == c), c("Amp","Amplen")])[,2])
        ## x <- n.var * sum(catlen)/totallen
        x <- (unique(vcf[which(as.character(vcf$Cat) == c), c("Amp","AmplenRatio")])[,2])
        x <- sum(x*n.var)
        obs <- which(as.character(vcf$Cat) == c)
        expcat <- rbind(expcat, c(c, x, "Expected"))
        expcat <- rbind(expcat, c(c, length(obs), "Observed"))
    }
    expcat <- data.frame(expcat)
    colnames(expcat) <- c("Category", "Counts", "Type")
    expcat$Counts <- as.numeric(as.character(expcat$Counts))
    expcat
}


############################################################
## Explore
############################################################

pdf("p20set.enrichment.pdf")
xmax <- max(vcfout$P20SizeInd)
## All
gg <- ggplot(vcfout, aes(x=P20SizeInd, fill=Cat)) +
    scale_x_discrete("Number of Hits", breaks=seq(1,xmax,1), labels = seq(1,xmax,1)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("All n=", length(unique(vcfout$Pos)),sep=""))
print(gg)
# enrichment
expcat <- countenrich(vcf)
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("All n=", length(unique(vcfout$Pos)),sep=""))
print(gg)
## Cov > 5000
fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 )
gg<-ggplot(vcfout[fltidx,], aes(x=P20SizeInd, fill=Cat)) +
    scale_x_discrete("Number of Hits", breaks=seq(1,xmax,1), labels = seq(1,xmax,1)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("MinCov > 5000", " n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
# enrichment
expcat <- countenrich(vcf[fltidx, ])
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("MinCov > 5000", " n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
dev.off()

## P20SizeInd >2
pdf("p20SizeInd.rep.pdf")
fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & vcf$Rep == "Rep")
gg<-ggplot(vcfout[fltidx,], aes(x=P20SizeInd, fill=Cat)) +
    scale_x_discrete("Number of Hits", breaks=seq(1,xmax,1), labels = seq(1,xmax,1)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("MinCov>5000 & Rep"," n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
expcat <- countenrich(vcf[fltidx,])
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("MinCov>5000 & Rep"," n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
dev.off()

## vcf$P20SizeInd >=2 & rank < 300
pdf("prank300.p20SizeInd.min2.pdf")
fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & vcf$P20SizeInd >=2 & vcf$PRank < 1000)
gg <- ggplot(vcfout[which(vcfout$P20SizeInd>=2),], aes(x=P20SizeInd, fill=Cat)) +
    scale_x_discrete("Number of Hits", breaks=seq(1,xmax,1), labels = seq(1,xmax,1)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Size>2 & Rank<1000"," n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
expcat <- countenrich(vcf[fltidx, ])
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("Size>2 & rank<1000"," n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
dev.off()


## vcf$P20SizeInd >=2
pdf("p20SizeInd.min2.pdf")
fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & vcf$P20SizeInd >=2)
gg <- ggplot(vcfout[which(vcfout$P20SizeInd>=2),], aes(x=P20SizeInd, fill=Cat)) +
    scale_x_discrete("Number of Hits", breaks=seq(1,xmax,1), labels = seq(1,xmax,1)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Hits identified in more than 1 sample"," n=", length(unique(vcfout$Pos[which(vcfout$P20SizeInd>=2)])), sep=""))
print(gg)
expcat <- countenrich(vcf[fltidx, ])
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("Hits identified in more than 1 sample"," n=", length(unique(vcfout$Pos[which(vcfout$P20SizeInd>=2)])), sep=""))
print(gg)
dev.off()

## top vs all
pdf("top.vs.all.pdf")
expcat <- countenrich(vcf)
gg <- ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("P20"," n=", length(unique(vcfout$Pos)), sep=""))
print(gg)
expcat.raw <- countenrich(anno_raw)
gg <- ggplot(expcat.raw, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("raw"," n=", length(unique(anno_raw$Pos)), sep=""))
print(gg)
dev.off()

#######################################################
#######################################################

badlibidx <- which(vcf$P20LibMax>=3 & vcf$P20SizeLib == 1)
topidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & vcf$RelPos>2 & (vcf$Rep == "Rep" | (vcf$PvlRank < 1000 & (vcf$P20SizeInd >=2))))
topidx <- topidx[!topidx %in% badlibidx]

## numbers
cat("bad library",length(badlibidx), "\n")
cat("coverage < 5000",length(which(vcf$MdCov<5000)), "\n")
cat("RelPos <= 2",length(which(vcf$RelPos<=2)), "\n")
cat("Rep",length(which(vcf$Rep == "Rep")), "\n")
cat("P20Size",length(which(vcf$P20SizeInd>=2)), "\n")

## dupsidx <- which(vcfout$genomicSuperDups != ".")
## topidx <- topidx[!topidx %in% dupsidx]

#######################################################
#######################################################


## relative position of these hits
pdf("posreltoedge.pdf")
ggplot(vcfout, aes(x=RelPos, fill=Cat)) +
    geom_histogram(aes(y=..count../sum(..count..))) +
    scale_x_discrete("Max distance to edge", breaks=seq(0,max(vcfout$RelPos[fltidx]),5), labels = seq(0,max(vcfout$RelPos[fltidx]),5)) +
            ggtitle("Variant position relative to Amplicon edge")
dev.off()

pdf("rate.top.pdf")
gg<-ggplot(vcfout, aes(x=log10(MxRate), fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle("All")
print(gg)
gg<-ggplot(vcfout[topidx,], aes(x=log10(MxRate), fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("MinCov>5000 & (PvlRank<100 or RateRank<100 & P20SizeInd>1)"," n=", length(unique(vcfout$Pos[fltidx])), sep=""))
print(gg)
dev.off()

vcfout$LJB23_SIFT_score_converted <- as.numeric(as.character(vcfout$LJB23_SIFT_score_converted))
vcfout$LJB23_SIFT_score_converted[which(is.na(vcfout$LJB23_SIFT_score_converted))] <- -1
vcfout$LJB23_Polyphen2_HVAR_score <- as.numeric(as.character(vcfout$LJB23_Polyphen2_HVAR_score))
vcfout$LJB23_Polyphen2_HVAR_score[which(is.na(vcfout$LJB23_Polyphen2_HVAR_score))] <- -1

anno_raw$LJB23_SIFT_score_converted <- as.numeric(as.character(anno_raw$LJB23_SIFT_score_converted))
anno_raw$LJB23_SIFT_score_converted[which(is.na(anno_raw$LJB23_SIFT_score_converted))] <- -1

anno_raw$LJB23_Polyphen2_HVAR_score <- as.numeric(as.character(anno_raw$LJB23_Polyphen2_HVAR_score))
anno_raw$LJB23_Polyphen2_HVAR_score[which(is.na(anno_raw$LJB23_Polyphen2_HVAR_score))] <- -1

pdf("sift.top.pdf")
gg<-ggplot(vcfout, aes(x=LJB23_SIFT_score_converted, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle("All")
print(gg)
gg<-ggplot(vcfout[topidx,], aes(x=LJB23_SIFT_score_converted, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Size>2 & Rank<300"," n=", length(unique(vcfout$Pos[topidx])), sep=""))
print(gg)
gg <- ggplot(anno_raw, aes(x=LJB23_SIFT_score_converted, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Input"," n=", dim(unique(anno_raw[,c("Chrm", "Pos")]))[1], sep=""))
print(gg)
dev.off()

pdf("polyphen.top.pdf")
gg<-ggplot(vcfout, aes(x=LJB23_Polyphen2_HVAR_score, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle("All")
print(gg)
gg<-ggplot(vcfout[topidx,], aes(x=LJB23_Polyphen2_HVAR_score, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Size>2 & Rank<300"," n=", length(unique(vcfout$Pos[topidx])), sep=""))
print(gg)
gg <- ggplot(anno_raw, aes(x=LJB23_Polyphen2_HVAR_score, fill=Cat)) +
        geom_histogram(aes(y=..count../sum(..count..))) +
            ggtitle(paste("Input"," n=", dim(unique(anno_raw[,c("Chrm", "Pos")]))[1], sep=""))
print(gg)
dev.off()

############################################################
## Output
############################################################
p20 <- sapply(strsplit(as.character(vcf$P20), ","), function(x){
    y <- as.numeric(x)
    y <- paste(as.character(y[y!=0]), collapse=";")
    y
})

p50 <- sapply(strsplit(as.character(vcf$P50), ","), function(x){
    if(length(x)>1){
        y <- as.numeric(x)
        y <- y[y!=0 & y!=1]
        y <- paste(as.character(y), collapse=";")
        y
    }else{
        "."
    }
})

vcfout$P20 <- p20
vcfout$P50 <- p50


write.table(vcfout, "anno.rank.txt", row.names=FALSE, sep="\t", quote=FALSE)
write.table(vcfout[topidx,], "anno.rank.top.txt", row.names=FALSE, sep="\t", quote=FALSE)

############################################################
### Singletons
############################################################

singletons <- read.csv("singleton.csv")

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
keys <- apply(singletons, 1, function(x){paste(trim(x),collapse="_")})
keysvcf <- apply(vcfout, 1, function(x){paste(trim(x[c(1,2,4,5,52)]),collapse="_")})

sidx <- which(keysvcf %in% keys)
write.table(vcfout[sidx,], "anno.rank.singletons.txt", row.names=FALSE, sep="\t", quote=FALSE)

                                        # write out sites
#write.table(unique(vcfout[topidx,c("Chr", "Pos")]), "/Users/zd1/volumn/wt/raindance/meta/topsites", row.names=FALSE, sep="\t", quote=FALSE)

############################################################
## Queries
############################################################

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
key="chr4_1803571_C_G"
keysvcf <- apply(vcfout, 1, function(x){paste(trim(x[c(1,2,4,5)]),collapse="_")})

match(key, keysvcf)

key2="chr2_25457243_G_A"
match(key2, keysvcf)


key3="chr11_533872_C_G"




############################################################
## Overlapping amplicons
############################################################

## # get overlapping and non overlapping amplicons
## duppos <- vcf[repi, c("Pos")]
## amppair <- unique(lapply(1:length(duppos), function(i){unique(vcf$Amp[which(vcf$Pos == duppos[i])])}))
## allamp <- as.character(unique(vcf$Amp))
## ampsingle <- allamp[! allamp %in% as.character(unlist(amppair))]

## prioritize ends here
