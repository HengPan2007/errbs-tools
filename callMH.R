#### Data input
sample <- '5485_05'
bamFile <- paste0('/home/hep2007/elementolab/Inghirami/data/ERRBS/bam/', sample, '.sorted.bam')
CGI <- read.table('/home/hep2007/lab/index/cpgIslandExt_hg19.txt')
CpG <- read.table(paste0('/home/hep2007/elementolab/Inghirami/data/ERRBS/CpG/cpg.', sample, '.mincov10.txt'), header=TRUE)
output <- paste0('/home/hep2007/lab/result/epiCGI/DLBCL/pattern.', sample, '.rds')
loci.num <- 4
loci.range <- 9
lowest.rc <- 60

#### load package ####
library(GenomicRanges)
library(Rsamtools)

#### Function ####
## Construct Methylation Pattern
patternList <- function(x){
  myPaste <- function(x){
    return(paste0(x[2], x[1]))
  }
  pattern <- c("C", "T")
  if(x==1) return(data.frame(pattern=pattern))
  else{
    while(x>1){
      pattern <- merge(c("C", "T"), pattern)
      pattern <- apply(pattern, 1, myPaste)
      x <- x -1
    }
    return(data.frame(pattern=pattern))
  }
}

##
patternMake <- function(lociPos, metPattern, bamFile, strand){
  pos <- as.numeric(as.data.frame(mcols(lociPos)))
  if (strand == 'F') {flag <- scanBamFlag(isMinusStrand=FALSE)} else if (strand == 'R') {flag <- scanBamFlag(isMinusStrand=TRUE)}
  param <- ScanBamParam(flag =flag, which=lociPos, what=c("pos", "seq"))
  tmp <- scanBam(bamFile, param=param)
  reads <- data.frame(Pos=tmp[[1]]$pos)
  reads$Seq <- as.character(tmp[[1]]$seq)
  index <- (reads$Pos <= pos[1]) & ((reads$Pos + nchar(reads$Seq) -1) >= pos[length(pos)])
  reads <- reads[index, ]
  if(length(reads$Pos)==0){return(rep(0, 16))}else{
    patternCollect <- function(x, strand){
      pattern <- DNAString(x[2])[pos-as.numeric(x[1])+1]
      if(strand == 'R'){pattern <- complement(pattern)}
      return(as.character(pattern))
    }
    reads$Pattern <- apply(as.matrix(reads), 1, patternCollect, strand)
    reads <- as.data.frame(table(reads$Pattern))
    colnames(reads)[1] <- "pattern"
    reads <- merge(metPattern, merge(metPattern, reads), all=TRUE)
    reads$Freq[is.na(reads$Freq)] <- 0
    return(reads$Freq)
  }
}

lociFilter <- function(anno){
  if (length(anno) == length(reduce(anno))){return(anno)} else{
    anno <- anno[order(mcols(anno)$rcs, decreasing=TRUE), ]
    tmp <- anno[2:length(anno)][!countOverlaps(anno[2:length(anno)], anno[1])]
    anno <- c(anno[1], lociFilter(tmp))
    return(anno)
  }
}
metPatternCGIsland <- function(idx, data, loci.num, lowest.rc){
  chr <- as.character(unique(data[idx, ]$chr))
  strand <- as.character(unique(data[idx, ]$strand))
  pos <- sort(data[idx, ]$base)
  pos.matrix <- cbind(pos[1:(length(pos)+1-loci.num)], pos[2:(length(pos)+2-loci.num)])
  if (loci.num > 2) {for (i in 3:loci.num){pos.matrix <- cbind(pos.matrix, pos[i:(length(pos)+i-loci.num)])}}
  pos.matrix <- as.data.frame(pos.matrix)
  pos.matrix <- pos.matrix[(pos.matrix[,4]-pos.matrix[,1]) <= loci.range, ]
  
  if (nrow(pos.matrix) == 0) {return(rep(0, 16))} else{
    loci.anno <- GRanges(seqnames=Rle(chr, nrow(pos.matrix)), ranges=IRanges(start=pos.matrix[,1], end=pos.matrix[,loci.num]), pos=pos.matrix)
    names(loci.anno) <- paste(chr, paste(start(loci.anno), end(loci.anno), sep='-'), sep=':')
    loci.anno.list <- split(loci.anno, names(loci.anno))
    result <- do.call(rbind, lapply(loci.anno.list, patternMake, metPattern, bamFile, strand))
    re.anno <- GRanges(seqnames=Rle(chr, nrow(pos.matrix)), ranges=IRanges(start=pos.matrix[,1], end=pos.matrix[,loci.num]), pattern=result, rcs=apply(result, 1, sum))
    re.anno <- re.anno[mcols(re.anno)$rcs > lowest.rc,]
    if(length(re.anno) == 0){return(rep(0, 16))} else{
      re.anno <- lociFilter(re.anno)
      result <- as.data.frame(mcols(re.anno)[1:ncol(mcols(re.anno))-1])
      rownames(result) <- paste(chr, paste(start(re.anno), end(re.anno), sep='-'), sep=':')
      return(result)
    }
  }
}
#### data passing ####
metPattern <- patternList(loci.num)
chr.list <- paste0("chr", 1:22)
chr.list[23] <- "chrX"
chr.list[24] <- "chrY"
chr.list[25] <- "chrM"

CpG.F <- CpG[CpG$strand=='F',]
CpG.R <- CpG[CpG$strand=='R',]
rm(CpG)

CGI.anno <- GRanges(seqnames=Rle(CGI$V2), ranges=IRanges(start=CGI$V3, end=CGI$V4))
seqlevels(CGI.anno, force=TRUE) <- chr.list
CpG.F.anno <- GRanges(seqnames=Rle(CpG.F$chr), ranges=IRanges(start=CpG.F$base, end=CpG.F$base))
CpG.R.anno <- GRanges(seqnames=Rle(CpG.R$chr), ranges=IRanges(start=CpG.R$base, end=CpG.R$base))

tmp <- countOverlaps(CGI.anno, CpG.F.anno)
tmp <- tmp >= loci.num
CGI.anno.F <- CGI.anno[tmp, ]
CGI.F.List <- as.list(findOverlaps(CGI.anno.F, CpG.F.anno))
tmp <- countOverlaps(CGI.anno, CpG.R.anno)
tmp <- tmp >= loci.num
CGI.anno.R <- CGI.anno[tmp, ]
CGI.R.List <- as.list(findOverlaps(CGI.anno.R, CpG.R.anno))

tmp <- do.call(rbind, lapply(CGI.F.List, metPatternCGIsland, data=CpG.F, loci.num=loci.num, lowest.rc=lowest.rc))
tmp <- tmp[apply(tmp, 1, sum) > lowest.rc, ]
colnames(tmp) <- metPattern$pattern
if (length(tmp) > 0) { rownames(tmp) <- paste('F', rownames(tmp), sep=':') }
result <- tmp
tmp <- do.call(rbind, lapply(CGI.R.List, metPatternCGIsland, data=CpG.R, loci.num=loci.num, lowest.rc=lowest.rc))
tmp <- tmp[apply(tmp, 1, sum) > lowest.rc, ]
colnames(tmp) <- metPattern$pattern
if (length(tmp) > 0) { rownames(tmp) <- paste('R', rownames(tmp), sep=':') }
result <- rbind(result, tmp)
saveRDS(result, file=output)