#########################################################
## export list
#########################################################

export_variants_within_2bp_of_primer <- function(vcf, badlibidx){
    iselect <- which(
    vcf$MdCov>5000 &
    vcf$RelPos<=2 &
    vcf$MxRate<0.03 & (vcf$Rep == "Rep" | vcf$P20SizeInd >=2))
    iselect <- setdiff(iselect, badlibidx)
    cat("length of selection:", length(iselect), "\n")

    vcf.1bp <- move_sample_cols_to_last(vcf[iselect, ])
    vcf.1bp <- convert_delimiter_to_semicolon(vcf.1bp)

    timestamp <- format(Sys.time(), "%Y_%b_%d_%H-%M-%S")

    write.table(
        vcf.1bp,
        file=paste0("2bp.otherwise_top_sites.", nrow(vcf.1bp), "list.", timestamp, ".tsv"),
        sep="\t", col.names=TRUE, row.names=FALSE)
}

export_toplist <- function(vcf){
    badlibidx <- which(vcf$P20LibMax>=3 & vcf$P20SizeLib == 1)
    badlibidx <- setdiff(badlibidx, which(vcf$Rep == "Rep"))

    iselect <- which(
        vcf$MdCov>5000 &
        vcf$RelPos>2 &
        vcf$MxRate<0.03 & (vcf$Rep == "Rep" | vcf$P20SizeInd >=2))
    iselect <- setdiff(iselect, badlibidx)
    cat("length of selection:", length(iselect), "\n")

    vcf.redo <- move_sample_cols_to_last(vcf[iselect, ])
    vcf.redo <- convert_delimiter_to_semicolon(vcf.redo)

    ## just check if the sites of interest are in the list
    stopifnot(
        nrow(vcf.redo[which(vcf.redo$Pos == 112888199 &
                            vcf.redo$Alt == "A"), metric_cols(vcf.redo)]) == 2)

    stopifnot(check_no_singletons(vcf.redo))
    stopifnot(check_know_replicates(vcf.redo))


    timestamp <- format(Sys.time(), "%Y_%b_%d_%H-%M-%S")

    write.table(
        vcf.redo,
        file=paste0("topsites.", nrow(vcf.redo), "list.", timestamp, ".tsv"),
        sep="\t", col.names=TRUE, row.names=FALSE)

}

#########################################################
## Hannah's toplist
#########################################################

get_hannahs_toplist <- function(){
    toplist.hr <- read.csv("RDanno.priority_24March15_Oct15_edit.csv",
                           stringsAsFactors=FALSE, sep="\t")
    toplist.hr$key <- get_key(toplist.hr)
    pcat("size of Hannah's toplist ", dim(toplist.hr), "\n")
    toplist.hr
}


compare_hannahs_toplist_against_local <- function(toplist.hr, vcf){
    toplist.hr.ft <-
    toplist.hr[which(!toplist.hr$key %in% badlib_labelled_within_the_463list$key),]

    ## look at the properties of the remaining sites
    toplist.hr.ft.vcf <- vcf[which(vcf$key %in% toplist.hr.ft$key),]

    ## pull out the singletons
    singletons <- toplist.hr.ft.vcf[which(toplist.hr.ft.vcf$P20SizeInd == 1),]
    write.table(
        singletons,
        file=paste0("singletons_", nrow(singletons), "list.tsv"),
        sep="\t", col.names=TRUE, row.names=FALSE)

    ## just check if data of single match between this list and the vcf obj
    vcf$P20SizeInd[which(vcf$key %in% toplist.hr.ft.vcf$key[which(toplist.hr.ft.vcf$P20SizeInd == 1)])]

    toplist.hr.ft.vcf.non_single <- toplist.hr.ft.vcf[which(!toplist.hr.ft.vcf$key %in% singletons$key),]

    setdiff(toplist.hr.ft.vcf$key,
            vcf[which(
                vcf$RelPos>2 &
                vcf$MdCov>5000 &
                vcf$MxRate<0.03 &
                vcf$RateRank < 1900 &
                vcf$PvlRank < 1300), "key"])

    cat("size of query (toplist.hr.ft.vcf.non_single)",
        nrow(toplist.hr.ft.vcf.non_single), "\n")

}

check_whats_dropped_from_hannahs_list <- function(){
    ########################################
    ## what's dropped from hannah's list
    ########################################
    ## also write out which ones will be dropped from the current list
    toplist.hr.ft.vcf.leftout <-
        toplist.hr.ft.vcf[
            which(toplist.hr.ft.vcf$key %in% setdiff(toplist.hr.ft.vcf$key, vcf.redo$key)), ]

    toplist.hr.ft.vcf.leftout.out <- move_sample_cols_to_last(toplist.hr.ft.vcf.leftout)

    write.table(
        toplist.hr.ft.vcf.leftout.out,
        file=paste0("sitesdropped_from_RDanno.priority_24March15_Oct15_edit.csv.",
                    nrow(toplist.hr.ft.vcf.leftout.out), "list.", timestamp, ".tsv"),
        sep="\t", col.names=TRUE, row.names=FALSE)

    cat("intersection:",
        length(intersect(toplist.hr.ft.vcf$key, vcf[iselect, "key"])),
        "\n")
}

####################################################
## check known stuff
####################################################

check_no_singletons <- function(invcf){
    known.singletons <- data.frame(do.call("rbind", strsplit(c("chr11_113934509_C", "chr11_113934514_C", "chr11_113934525_C", "chr11_113934547_A", "chr11_113934548_G", "chr11_113934563_C", "chr11_113934581_C", "chr11_113934606_G", "chr11_113934631_T", "chr11_113934640_C", "chr11_113934667_C", "chr11_113934689_C", "chr11_113934718_A", "chr11_113934734_C", "chr11_113934739_G", "chr11_113934697_T", "chr11_113934706_C", "chr11_113934728_G", "chr11_113934732_C", "chr11_113934738_A", "chr11_113934752_C", "chr11_113934758_C", "chr11_113934760_T", "chr11_113934762_C", "chr11_113934776_T", "chr11_113934777_A", "chr11_113934782_A", "chr11_113934786_G", "chr11_113934801_G", "chr11_113934806_A", "chr11_113934807_A"), split="_")), stringsAsFactors=FALSE)
    colnames(known.singletons) <- c("Chrm", "Pos", "Alt")
    found <- merge(known.singletons,
                   invcf,
                   by.x=c("Chrm", "Pos", "Alt"),
                   by.y=c("Chrm", "Pos", "Alt"))
    nrow(found) == 0
}


check_know_replicates <- function(invcf){
    known.replicates <- data.frame(
        do.call(
            "rbind",
            strsplit(c("chr12_112940006_A", "chr12_112888198_C"), split="_")
        ),
        stringsAsFactors=FALSE
    )
    colnames(known.replicates) <- c("Chrm", "Pos", "Alt")
    nfound <- nrow(merge(known.replicates,
                         invcf,
                         by.x=c("Chrm", "Pos", "Alt"),
                         by.y=c("Chrm", "Pos", "Alt")))
    nfound > 0
}


get_replicated_sites <- function(vcf){
    ## Replicated - a variant identified in overlapping amplicons in *same* samples
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

    known.singletons <- data.frame(do.call("rbind", strsplit(c("chr11_113934509_C", "chr11_113934514_C", "chr11_113934525_C", "chr11_113934547_A", "chr11_113934548_G", "chr11_113934563_C", "chr11_113934581_C", "chr11_113934606_G", "chr11_113934631_T", "chr11_113934640_C", "chr11_113934667_C", "chr11_113934689_C", "chr11_113934718_A", "chr11_113934734_C", "chr11_113934739_G", "chr11_113934697_T", "chr11_113934706_C", "chr11_113934728_G", "chr11_113934732_C", "chr11_113934738_A", "chr11_113934752_C", "chr11_113934758_C", "chr11_113934760_T", "chr11_113934762_C", "chr11_113934776_T", "chr11_113934777_A", "chr11_113934782_A", "chr11_113934786_G", "chr11_113934801_G", "chr11_113934806_A", "chr11_113934807_A"), split="_")), stringsAsFactors=FALSE)
    colnames(known.singletons) <- c("Chrm", "Pos", "Alt")

    vcf.rep.known.singletons <- merge(known.singletons,
                                      vcf.rep,
                                      by.x=c("Chrm", "Pos", "Alt"),
                                      by.y=c("Chrm", "Pos", "Alt"))

    stopifnot(nrow(vcf.rep.known.singletons) == 0)
    vcf.rep

}

####################################################
## helper functions
####################################################

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

## back up the original dat
## vcf.bk <- vcf
## vcfout.bk <- vcfout

metric_cols <- function(dat){
    colnamestoskip <- c(
	1:288,
	paste("X", 1:288, sep=""),
	colnames(anno)
    )
    unique(c(which(!colnames(dat) %in% colnamestoskip),
	     which(colnames(dat) %in% c("Ref", "Alt", "Rep"))))
}

move_sample_cols_to_last <- function(invcf){
    innames <- colnames(invcf)
    nnames <- length(innames)
    start <- which(innames == "1")
    end <- which(innames == "288")
    if(end == nnames) {
        ordered_cols <- c(innames[1:(start-1)], start:end)
    } else {
        ordered_cols <- c(innames[1:(start-1)], innames[(end+1):nnames], 1:288)
    }
    invcf[, ordered_cols]
}


convert_delimiter_to_semicolon <- function(invcf){
    invcf$P20 <- sapply(invcf$P20, function(x){gsub(",", ";", x)})
    invcf$P50 <- sapply(invcf$P50, function(x){gsub(",", ";", x)})
    invcf
}
