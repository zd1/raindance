
library(fpc)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(gplots)
library(GMD)
library(plyr)

#####################
# when run in batch
#####################
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
ampkey <- args[1]
wkdir <- args[2]
sample.anno <- read.csv("/users/mcvean/zd1/volumn/raindance/meta/sample_seq.csv", head=TRUE)
rm(args)

##################
### when testing
##################
## ## ampkey <-  "50:FGFR2_5:chr10:123279598-123279714"
## ## ampkey <-  "48:FGFR2_5:chr10:123279481-123279643"
## ## ampkey <- "49:FGFR2_5:chr10:123279547-123279684"
## ampkey <-  "174:KRAS:chr12:25398189-25398313"
## #ampkey <-  "175:KRAS:chr12:25398214-25398362"
#wkdir <- paste("/Users/zd1/cloud/data/raindance/pileup/sum/", ampkey, sep="")

## ampkey <- "229:TP53_11:chr17:7579573-7579729"
## wkdir <- paste("/Users/zd1/volumn/wt/raindance/sum/", ampkey, sep="")
## sample.anno <- read.csv("/Users/zd1/cloud/data/raindance/meta/sample_seq.csv", head=TRUE)


######################
## Load data
######################

setwd(wkdir)
if (!file.exists(wkdir)){
    warning("Directory doesn't exist, amp proobably failed")
    quit()
}

ampfile <- paste(ampkey, ".csv",sep="")
if (!file.exists(ampfile)){
    warning("Amplicon csv doesn't exist, amp proobably failed")
    quit()
}

data <- read.csv(ampfile, head=TRUE, colClasses=c("character",  "character", "numeric", "character",  rep("numeric", 288*4)))

# get annotation, match it with samples read from data
seqsamplelist <- colnames(data)[-c(1,2,3,4)][seq(1,1440,5)]
samples <- sapply(1:length(seqsamplelist), function(i){paste(strsplit(seqsamplelist[i], "_")[[1]][1:3], collapse="_")})

allsamples <- sample.anno[match(samples, sample.anno$seqsample),]

ref <- data[,4]
pos <- data[,3]
chrm <- unique(data[,2])
n.sample <- nrow(allsamples)
n.meta <- 4 # meta columns
n.site<-nrow(data)

libcol <- allsamples$Lib
indcol <- allsamples$individual
allind <- unique(indcol)
n.ind <- length(allind)
bccol <- allsamples$BCNumber

# store index for the highest and second highest allele
consensus.o <- array(0, c(n.site, 2))
ncr <- array(0, n.site,)
highestsecond <- array(0, n.site)
ncc.sample <- array(0, c(n.site, n.sample))
ncr.sample <- array(0, c(n.site, n.sample))

# separate counts for ATGC
allele.sample <- array(0, c(n.site, n.sample, 4))
allele.sample.c <- array(0, c(n.site, n.sample, 4))

a.let<-c("A", "C", "G", "T");
snp.thresh<-0.1;
n.snps <-  0
snp.sites <- c()
site.qual <- array(0,c(n.site,n.sample))

for (s in 1:n.site) {
    sub <- data[s, c(-1,-2,-3,-4)]
    d.s <- array(as.numeric(data[s, 5:ncol(data)]), c(5, n.sample))
    site.qual[s, ] <- d.s[5,]
    d.s <- d.s[1:4,]
    ts<-apply(d.s, 2, sum, na.rm=T); # sum across rows to get total counts for each sample
                                        # we found one sample that has got its ref count set to 0, use average for now
                                        # to be removed
    ts[which(ts<1000)] = mean(ts)
    o <- apply(d.s,2,function(x){order(x, decreasing=T)})
    fs <- aaply(d.s, 1, "/", ts) # get allele frequency
    l.highest <- round(mean(o[1,]))
    l.second <- round(mean(o[2,]))
    l.highchunk <- which.max(sapply(1:dim(o)[2], function(x){fs[o[2,x],x]}))
                                        # the most elevated allele
    if (mean(fs[l.second,])> snp.thresh){
        consensus.o[s, 2]<- l.highest;
    }
    ncc.sample[s,] <-  apply(d.s[c(-l.highest),],2, sum)
    ncr.sample[s,] <-  apply(fs[c(-l.highest),],2, sum)
    # store counts for each base separately, set ref = 0
    allele.sample[s,,] = t(fs)
    allele.sample.c[s,,] = t(d.s)
    allele.sample[s,, which(a.let == ref[s])] = 0
    allele.sample.c[s,, which(a.let == ref[s])] = 0
    if (consensus.o[s,2]==0) {
        count.consen = signif(100*mean(apply(fs[-l.highest,],2,sum)), 3)
        cat("\nSite ", s, " | Hom ", a.let[l.highest], "\t| AverageNonConRead = ", count.consen ,  "%", sep="");
        ncr[s]<-median(apply(fs[-l.highest,],2,sum));
    } else { # if it's a het
        count.het.combined = signif(100*mean(apply(fs[c(-l.highest, -l.second),],2,sum)), 3)
        cat("\nSite ", s, " | Het ", a.let[l.highest], "\t| AverageNonConRead = ", count.het.combined,  "%", sep="");
        ncr[s]<-median(apply(fs[c(-l.highest, -l.second),],2,sum));
        n.snps<-n.snps+1;
        snp.sites <- c(snp.sites, s)
        # we ignore SNPs
        allele.sample[s,,] = 0
        allele.sample.c[s,,] = 0
    }
}


######################
## Testings
######################

# pca then gaussian regression
removepc <- function(samplebysite, pc){
    fit <- princomp(samplebysite, cor=FALSE)
    gen.10pc <- fit$scores[,1:pc] %*% t(fit$loadings[,1:pc])
    sample.pc <- samplebysite-gen.10pc
    sample.pc
}

# enforce upper limit. heatmap.3 doesn't have an upper limit for value range
forceRange <- function(dataset, maxv = 0.003){
    dataset[which(dataset > maxv)] = maxv
    dataset
}

# regress out known covariates
glmres <- function (samplebysite){
    res <- array(0, dim=dim(samplebysite))
    for(i in 1:(dim(samplebysite)[2])){
        qual <- site.qual[i,]
        qual[which(is.na(qual))] <- median(qual, na.rm=TRUE)
        qual[!is.finite(qual)] <- median(qual[which(is.finite(qual) & !is.na(qual))])
        glmfit <- glm(samplebysite[,i] ~ qual + allsamples$"Flow.cell" + allsamples$Lib + allsamples$individual, family=gaussian())
        res[,i] <- residuals(glmfit)
    }
    res
    ## gs <- sapply(1:dim(samplebysite)[2], function(i){
    ##     glmfit <- glm(samplebysite[,i] ~ site.qual[i,] + allsamples$"Flow.cell" + allsamples$Lib + allsamples$individual, family=gaussian())
    ##     residuals(glmfit)
    ## })
}

# variance stablizer by taking the sd of the values below the median
vstb <- function(x){
    sqrt(mean((x[which(x<median(x))]-median(x))^2))
}


# remove median of each lane
laneKnn <- function(samplebysite){
    res <- array(0, dim=dim(samplebysite))
    n.site <- dim(samplebysite)[2]
    lanes <- unique(allsamples$Illumina.lane)
    for(i in 1:n.site){
        lnmed <- c()
        lnsd <- c()
        for(l in lanes){
            li <- which(allsamples$Illumina.lane == l)
            lnmed <- c(lnmed, median(samplebysite[li,i]))
            lnsd <- c(lnsd, vstb(samplebysite[li,i]))
        }
        lnsdscaled <- lnsd/mean(lnsd)
        for(ll in 1:length(lanes)){
            li <- which(allsamples$Illumina.lane == lanes[ll])
            res[li,i] <- (samplebysite[li,i] - lnmed[ll])
        }
    }
     res
}

# normalise by inter-quantile-range
interqt <- function(samplebysite){
    x <- samplebysite
    q1<-apply(x, 2, quantile, 0.25)
    q2<-apply(x, 2, quantile, 0.75)
    zz <- x
    for (i in 1:ncol(x)) zz[,i]<-x[,i]/(q2[i]-q1[i])
    zz
}

# priori. if the max count at a site is below 10 we don't believe it's mutation
ptest_sitebg <- function(samplebysite, allelei, mincount=10){
    ## sd0 <- sd(x, na.rm=T)
    ## mu0 <- median(x, na.rm=T)
    mut.pvl <- array(0, c(n.sample, n.site))
    for(s in 1:n.site){
        if(max(allele.sample.c[s,, allelei], rm.na=TRUE)<mincount){
            mut.pvl[,s] = 1
            next
        }
        refi <- which(a.let == ref[s])
        if(a.let[allelei] == ref[s]){
            mut.pvl[,s] = 1
        }else{
            q1 <- quantile(samplebysite[,s], probs=0.25)
            q3 <- quantile(samplebysite[,s], probs=0.75)
            mu0 <- median(samplebysite[,s])
            sd0 <- q3 - q1
            mut.pvl[,s] <- pnorm(samplebysite[,s], mean=mu0, sd=sd0, lower.tail=FALSE)
            #mut.pvl[,s] <- pnorm(samplebysite[,s], lower.tail=FALSE)
        }
    }
    mut.pvl
}


call_list <- function(pvl, allelei, mincount=10, bonf=TRUE){
    gg.mx <- cbind(melt(pvl), a.let[allelei])
    colnames(gg.mx) <- c("Sample", "Site", "P", "Allele")
    gg.mx$Exp <- rank(gg.mx$P)/(dim(gg.mx)[1]+1)
    gg.mx$bonf <- p.adjust(gg.mx$P, method = "bonferroni")
    gg.mx$fdr <- p.adjust(gg.mx$P, method = "fdr")
    if(bonf){
        threshold <- which(gg.mx$bonf < 0.01) # bonf at 1%
    }else{
        threshold <- which(gg.mx$fdr <0.01) # fdr at 1%
    }
    res <- gg.mx[threshold, ]
    if(dim(res)[1] == 0){
        return(NULL)
    }
    res$SamName <- allsamples$seqsample[res[,"Sample"]]
    res$Chrm <- chrm
    res$Pos <- pos[res[,"Site"]]
    res$Rate <- sapply(1:nrow(res), function(i){allele.sample[res[i,"Site"],res[i,"Sample"], allelei]})
    res$Count <- sapply(1:nrow(res), function(i){allele.sample.c[res[i,"Site"],res[i,"Sample"], allelei]})
    res$Ref <- ref[res$Site]
    res <- res[which(res$Count>mincount), ]
    if(dim(res)[1] == 0){
        return(NULL)
    }else{
        return(res)
    }
    res$NCC <- sapply(1:nrow(res), function(i){ncc.sample[res[i,"Site"],res[i,"Sample"]]})
    res$NCR <- sapply(1:nrow(res), function(i){ncr.sample[res[i,"Site"],res[i,"Sample"]]})
    if(dim(res)[1] == 0){
        return(NULL)
    }else{
        return(res)
    }
}

pval_plots <- function(allp){
    gg.mx <- c()
    for(a in 1:4){
        gg.mx <- rbind(gg.mx, cbind(melt(allp[,,a]), a.let[a]))
    }
    colnames(gg.mx) <- c("Sample", "Site", "P", "Allele")
    gg.mx$Exp <- rank(gg.mx$P)/(dim(gg.mx)[1]+1)
    gg.mx$P[which(gg.mx$P==0)] <- 1E-290 # ones with 0 pvalues
    g1 <- ggplot(gg.mx, aes(x=Site, y=-log10(P), colour=Allele)) +
          geom_point()  + ggtitle(ampkey)
    sigi <- which(gg.mx$P < 0.01/(dim(gg.mx)[1])*3) # three different types of mut
    g2 <- ggplot(gg.mx, aes(x=-log10(Exp), y=-log10(P), colour=Allele)) +
            geom_point() +
            ## annotate("text", x = -log10(gg.mx[sigi, "Exp"])+0.05, y = -log10(gg.mx[sigi,"P"]),
            ##          label = pos[gg.mx[sigi, "Site"]], size=3) +
                         geom_abline(intercept = 0, slope = 1, colour="green", size=1) 
                             #geom_hline(aes(-log10(0.05/dim(gg.mx)[1]), colour="Bonf"))
    grid.arrange(g1, g2, ncol=2)
}

natoval <- function(x, val=0){
    x[is.na(x)] <- val
    return(x)
}

#heatmap.3(t(allele.sample[,,1]), dendrogram = "none", Rowv=FALSE, Colv=FALSE, main=paste(ampkey,"raw"))

allp <- array(0, dim=c(n.sample, n.site, 4))
lambda <- array(0, dim=c(n.site, 4))
q60 <- array(0, dim=c(n.site, 4))
q75 <- array(0, dim=c(n.site, 4))

hits <- c()

for(a in a.let){
    allelei <- which(a.let == a)
    mincount <- 10
    xglm <- glmres(t(allele.sample[,, allelei]))
    xknn <- laneKnn(xglm)
    x <- interqt(xknn)
    ld <- natoval(apply(x,2,median))
    qq60 <- natoval(apply(x,2,function(x){quantile(x, probs=0.6, na.rm=TRUE)/qnorm(0.6)}))
    qq75 <- natoval(apply(x,2,function(x){quantile(x, probs=0.75, na.rm=TRUE)/qnorm(0.75)}))
    lambda[,allelei] <- ld
    q60[,allelei] <- qq60
    q75[,allelei] <- qq75
    pvl <- ptest_sitebg(x, allelei, mincount)
    ## sig <- call_list(pvl, allelei, bonf=FALSE)
    ## if(!is.null(sig)){
    ##     hits <- rbind(hits, sig)
    ## }
    allp[,, allelei] <- pvl
}

# plot p vals
pdf("p.glm.knn.ln.itq.qtsd.pdf", width=12, height=6)
pval_plots(allp)
dev.off()

###################################################
## Output everything
###################################################
save.image(paste(ampkey,".Rdata",sep=""))

for(a in 1:4){
    meta <- data.frame(Chrm=chrm, Pos=pos, Ref=ref,
                       Alt=rep(a.let[a],length(ref)),
                       AvgQual = round(apply(site.qual,1,mean), 3),
                       AmpPos = seq(length(pos)),
                       AmpPosNeg = seq(from=length(pos), to=1),
                       Lambda=signif(lambda[,a], 3), 
                       Q60=signif(q60[,a], 3),
                       Q75=signif(q75[,a], 3)
                       )
    count <- allele.sample.c[,,a]
    rate <-  allele.sample[,,a]
    colnames(rate) <- samples
    rc <- matrix(paste(round(t(-log10(allp[,,a])),5), count, round(rate,8), sep=";"), 
                 nrow=nrow(rate), dimnames=dimnames(rate))
    write.csv(cbind(meta, rc), file = paste(ampkey,".",a.let[a],".csv",sep=""),
          row.names = FALSE, quote=FALSE)
}

colnames(ncc.sample) <- samples
write.csv(cbind(Pos=pos, ncc.sample), file = paste(ampkey,".ncc.csv",sep=""),
          row.names = FALSE, quote=FALSE)
colnames(ncr.sample) <- samples
write.csv(cbind(Pos=pos, ncr.sample), file = paste(ampkey,".ncr.csv",sep=""),
          row.names = FALSE, quote=FALSE)


###################################################
## Write out putative hits
###################################################

## if(length(hits) > 0){
##     # organise output formats
##     hits$Amp <- ampkey
##     hits <- hits[, c("Amp","Chrm", "Pos", "SamName", "Sample", "Site",  "Ref",  "Allele",
##                      "P", "Exp", "fdr", "bonf", "Rate", "Count",  "NCC", "NCR")]
##     hits$Depth <- hits$NCC/hits$NCR
##     write.table(hits, paste(ampkey,".calls.qtsd.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
## }


###################################################
## debugging as some sample/site were not called
###################################################

## tpos <- 25398281
## knowsam <- c(274, 115,116,118,120,121,122,124, 126)

## ## tpos <- 25398284
## ## knowsam <- c(210, 218, 246, 90, 92, 93, 96, 98, 46, 60, 259) 

## ipos <- which(pos==tpos)
## calledsams <- hits$Sample[which(hits$Pos==tpos)]

## sample.v <- data.frame(Sample=allsamples$seqsample,
##                        Lane = allsamples$Illumina.lane,
##                        ## Lib = allsamples$Lib,
##                        ## FC = allsamples$"Flow.cell",
##                        ## Ind = allsamples$individual,
##                        ## Slide = allsamples$slide.ID,
##                        ORI=allele.sample[ipos,, 4],
##                        GLM=glmres(t(allele.sample[,, 4]))[,ipos],
##                        KNN=laneKnn(glmres(t(allele.sample[,, 4])))[,ipos],
##                        QT=interqt(laneKnn(glmres(t(allele.sample[,, 4]))))[,ipos],
##                        P=-log10(allp[,ipos,4]))

## sample.v$anno <- "null"
## sample.v$anno[knowsam] <- "hits"
## sample.v$anno[calledsams] <- "called"
## sample.v$anno <- as.factor(sample.v$anno)

## pdf(paste("pairs.remove.lane.effect.",tpos,".pdf", sep=""))
## pairs(sample.v[,-8], cex=0.7, pch=16, col=sample.v$anno, main="25398281")
## dev.off()




## sample.t.ori <- melt(t(allele.sample[,, 4]))
## sample.t.glm <- melt(glmres(t(allele.sample[,, 4])))
## sample.t.knn <- melt(laneKnn(glmres(t(allele.sample[,, 4]))))
## sample.t.qt <- melt(interqt(laneKnn(glmres(t(allele.sample[,, 4])))))

## sample.t <- cbind(sample.t.ori, sample.t.glm[,3], sample.t.knn[,3], sample.t.qt[,3])
## colnames(sample.t) <- c("Sample", "Site", "Rate", "GLM", "KNN", "QT")
## sample.t <- sample.t[which(sample.t$Site == hitpos),]
## sample.t$Count <- sapply(1:nrow(sample.t), function(i){allele.sample.c[sample.t[i,"Site"],sample.t[i,"Sample"], 4]})
## sample.t$logP <- sapply(1:nrow(sample.t), function(i){-log10(allp[sample.t[i,"Sample"],sample.t[i,"Site"], 4])})
## #anno.which <- which(sample.t$NormRate > 0.002)
## sample.t$hit <- "green"
## sample.t$hit[sample.t$Sample %in% calledsams] <- "blue"
## sample.t$hit[sample.t$Sample %in% knowsam] <- "red"
## sample.t$SamName <- allsamples$seqsample[sample.t[,"Sample"]]
## sample.t <- cbind(sample.t, allsamples[match(sample.t$SamName, allsamples$seqsample),])
## sample.t <- sample.t[, which(!colnames(sample.t) %in% c("X", "SamName", "BCNumber", "Lib_BC", "Testis.Sample", "Testis.slide.ID", "X.1", "Ind.Slide.sample", "seqsample", "BC",  "Site")) ]

## pdf("diagnosis.pair.25398281.pdf", width=8, height=6)
## pairs(sample.t, col=sample.t$hit, cex=0.7, pch=16)
## dev.off()


## ggplot(sample.t, aes(x=Sample, y=logP, colour=hit)) + geom_point()


    ## annotate("text", x=sample.t$Site[anno.which] + 2, y = sample.t$NormRate[anno.which], label = sample.t$Sample[anno.which], size=3)
