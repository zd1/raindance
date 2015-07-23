
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
sample.anno <- read.csv("meta/sample_seq.csv", head=TRUE)
rm(args)

######################
## Testing
######################

ampkey <-  "50:FGFR2_5:chr10:123279598-123279714"
wkdir <- "/Users/zd1/cloud/data/raindance/pileup/sum/50:FGFR2_5:chr10:123279598-123279714/pub"
sample.anno <- read.csv("/Users/zd1/cloud/data/raindance/meta/sample_seq.csv", head=TRUE)


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
    gg.mx$pos <- pos[gg.mx$Site]
    g1 <- ggplot(gg.mx, aes(x=pos, y=-log10(P), colour=Allele)) +
        geom_point()  + ggtitle(ampkey) + xlab("Position")
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
## heatmaps
###################################################

poslab <- rep("", length(pos))
poslab[seq(1,length(pos),10)] <- pos[seq(1,length(pos),10)]

lane_colors = rep(rep(c("darkorchid", "darkred"), each=16), 9)
fc_colors = rep(c("green","darkgreen"), each=144)

rlab=t(cbind(lane_colors, fc_colors))

rownames(rlab)=c("Lane","Flow Cell")

pdf("pileup.c.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(t(allele.sample[,,2]),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()



################################################
## coverage
################################################
cvg <- read.table("coverages.txt")
itarget <- read.table("target_index.txt")

postt <-  pos[itarget[,1]]
poslab <- rep("", length(postt))
poslab[seq(1,length(postt),10)] <- postt[seq(1,length(postt),10)]
lane_colors = rep(rep(c("darkorchid", "darkred"), each=16), 9)
fc_colors = rep(c("green","darkgreen"), each=144)
rlab=t(cbind(lane_colors, fc_colors))
rownames(rlab)=c("Lane","Flow Cell")


pdf("coverage.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(cvg[,itarget[,1]],
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()


################################################
## Normalisation stepwise
################################################

allelei <- which(a.let == "C")
ori <- t(allele.sample[,, allelei])
xglm <- glmres(ori)
xknn <- laneKnn(xglm)

pdf("norm.ori.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(forceRange(ori),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()


pdf("norm.glm.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(forceRange(xglm),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()

pdf("norm.knn.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(forceRange(xknn),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()



################################################
## NCC
################################################

pdf("ncc.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(t(ncc.sample),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Count")
dev.off()

pdf("qual.mx.pdf", width=10, height=8)
main_title=""
par(cex.main=1)
heatmap.3(t(site.qual),
          Rowv = FALSE,
          Colv = FALSE,
          scale="none",
          dendrogram="none",
          margins=c(6,12),
          RowSideColors=rlab,
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none",
          main=main_title,
          labCol=poslab,
          ## labRow=seq(1,288),
          col=topo.colors(200),
          labRow=FALSE,
          cexRow=1,
          cexCol=0.8,
          RowSideColorsSize=2,
          KeyValueName="Phred Q")
dev.off()





################################################
## customised heatmaps
################################################
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }
 
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
 
    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
 
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
 
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}



