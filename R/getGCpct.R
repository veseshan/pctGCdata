getGCpct <- function(chrom, pos, gbuild=c("hg19", "hg38", "hg18", "mm9", "mm10")) {
    gbuild <- match.arg(gbuild)
    # check that the chromosome is valid
    if (length(chrom) > 1) {
        warning("only the first element of chrom is used")
        chrom <- chrom[1]
    }
    if (gbuild %in% c("hg19", "hg38", "hg18")) {
        if (is.character(chrom)) {
            chrom <- match(chrom, c(1:22, "X", "Y", "MT"))
        } else {
            chrom <- match(chrom, 1:25)
        }
    } else {
        if (is.character(chrom)) {
            chrom <- match(chrom, c(1:19, "X", "Y", "MT"))
        } else {
            chrom <- match(chrom, 1:22)
        }
    }
    if (is.na(chrom)) stop("incorrect chrom sepcification")
    # if the needed genome build is in workspace set it to gcpctdb
    gcpctdb <- get(gbuild, pos="package:pctGCdata")
    # gc percentages of specified chromosome
    gcpct <- gcpctdb[[chrom]]
    # left and right intervals for the genomic positions
    jleft <- ceiling((pos-499.5)/100)
    jright <- ceiling((pos-399.5)/100)
    # make sure that indices are valid
    jleft[jleft < 1] <- 1
    jleft[jleft > length(gcpct)] <- length(gcpct)
    jright[jright < 1] <- 1
    jright[jright > length(gcpct)] <- length(gcpct)
    # difference between pos and midpoints
    jj <- 100 - (pos %% 100)
    # linear combination of the two GC percentages
    round((jj*gcpct[jleft] + (100-jj)*gcpct[jright])/100, digits=3)
}
