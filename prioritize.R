
## [[file:~/dev/org/rain.org::*prioritize%20variants][prioritize]]

#library(caret)
library(ggplot2)
library(reshape)

wkdir="/Users/zd1/cloud/data/raindance/pileup/sum/calls/annovar"
setwd(wkdir)
samples <- read.table("/Users/zd1/cloud/data/raindance/meta/sample_seq.csv", sep=",", header=TRUE)
samplelist <- read.table("/Users/zd1/cloud/data/raindance/meta/samplelist", sep=",", header=FALSE)

samples$id <- match(samples$seqsample, samplelist$V1)


vcf=read.table("annovar.vcf.gz.avinput", header=TRUE, sep="\t", comment.char="")
colnames(vcf)[1:14] <- c("Chrm", "Pos", "Pos", "Ref", "Alt",
                         "Amp", "MaxP", "P20", "P50", "Phred",
                         "AC_NFE", "AC", "AF", "CNT")

colnames(vcf)[15:302] <- seq(1,288)
anno=read.table("annovar.vcf.gz.annovar.out.hg19_multianno.no2ndrow.txt", header=TRUE, sep="\t")
ampcat=read.table("/Users/zd1/cloud/data/raindance/meta/ampcat.csv", header=FALSE, sep=",")
colnames(ampcat) <- c("AmpID", "Amp", "CatNum")



# annotation for the entire input region covered by this experiment
anno_raw=read.table("raw.vcf.gz.annovar.out.hg19_multianno.no2ndrow.txt", header=TRUE, sep="\t")
anno_avinput=read.table("raw.vcf.gz.avinput", header=TRUE, sep="\t", comment.char="")
colnames(anno_avinput) <- c("Chrm", "Pos", "Pos1", "Ref", "Alt", "Amp", "MaxP", "P20", "P50")
anno_raw <- cbind(anno_raw, anno_avinput)

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


############################################################################
## Mutation Pattern
############################################################################

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

amplen <- sapply(1:dim(anno_raw)[1], function(x){
    y <- as.integer(strsplit(strsplit(as.character(anno_raw$Amp)[x], ":")[[1]][4],"-")[[1]])
    abs(y[1] - y[2])
})
anno_raw$Amplen <- amplen


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

anno_raw$Cat <- ampcat$Cat[match(as.character(anno_raw$Amp), as.character(ampcat$Amp))]
anno_raw$AmplenRatio <- ampcat$AmplenRatio[match(as.character(anno_raw$Amp), as.character(ampcat$Amp))]

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
a.let <- c("A", "C", "G", "T")
altcol <- 5
sampstart <- 15
n.sam <- 288

vcfmeta <- apply(vcf, 1, function(x){
    vsi <- which(a.let == as.character(x[5])[[1]])
    x <- x[15:302] # samples
    vsc <- lapply(x,function(y){strsplit(strsplit(as.character(y),":")[[1]][1], ",")})
    mdcov <- median(sapply(vsc, function(z){sum(as.numeric(z[[1]]))}))
    mxrate <- max(sapply(vsc, function(z){as.numeric(z[[1]])[vsi]/sum(as.numeric(z[[1]]))}), na.rm=TRUE)
    c(mdcov,mxrate)
})

vcf$MdCov <- vcfmeta[1,]
vcf$MxRate <- vcfmeta[2,]
vcf$logRate <- log10(vcf$MxRate)
vcf$PvlRank <- ave(-vcf$MaxP, vcf$Cat, FUN=rank)
vcf$RateRank <- ave(-vcf$logRate, vcf$Cat, FUN=rank) 

vcf$P20Size <- sapply(strsplit(as.character(vcf$P20),","), length)
vcf$P20SizePvlRank <- ave(-vcf$MaxP, vcf$P20Size, FUN=rank) 

vcf$PRank <- rank(-vcf$MaxP)

# Replicated
repi <- which(duplicated(vcf[,c("Chrm", "Pos", "Alt")]) | duplicated(vcf[,c("Chrm", "Pos", "Alt")], fromLast = TRUE))
vcf$Rep <- "Uni"
vcf$Rep[repi] <- "Rep"
vcf$Rep <- as.factor(vcf$Rep)

# libraries
libs <- unique(samples$Illumina.lane)
libc <- array(0, c(dim(vcf)[1], length(libs)))

# individuals
inds <- unique(samples$individual)
indc <- array(0,c(dim(vcf)[1],length(inds)))

# add counts              
for(i in seq(1,dim(vcf)[1])){
    sam_i <- samples[match(as.integer(strsplit(as.character(vcf$P20[i]), ",")[[1]]), samples$id),]
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

vcf$P20SizeInd <- apply(indc,1,max)
vcf$P20SizeLib <- apply(libc, 1, function(x){length(x[which(x>0)])})
vcf$P20LibMax <- apply(libc,1,max)

features<-c("Chrm","Pos","Ref"
,"Alt","Amp"
,"MaxP","P20"
,"P50","Phred"
,"AC_NFE","AC"
,"AF","CNT","Cat"
,"MdCov","MxRate"
,"logRate","Rep", "RelPos"
,"RateRank","P20Size", "P20SizeInd", "P20SizeLib" ,"P20SizePvlRank", "PvlRank", "Amplen", "Mut", seq(1,288))

vcfout <- cbind(anno, indc, libc, vcf[,features])

#save.image("vcf.added.Rdata")

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
## Enrichment
############################################################

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

fltidx <- which(vcf$MdCov>5000 & vcf$MxRate>0 & (vcf$PvlRank<100 | vcf$RateRank<100 | vcf$Rep == "Rep") &  vcf$P20Size>1)
n <- length(unique(vcfout$Pos[fltidx]))
expcat <- countenrich(vcf[fltidx,])
pdf(paste("expcounts_top.",n,".pdf",sep=""))
ggplot(expcat, aes(x=Category, y=Counts,fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
        ggtitle(paste("MinCov>5000 & (PvlRank<100 or RateRank<100 or Rep) & Size>1"," n=",n, sep=""))
dev.off()

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
