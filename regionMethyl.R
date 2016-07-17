#### arguments preparation
args <- (commandArgs(TRUE))
input_dir <- ''
output_dir <- '.'
regions <- ''

if (length(args) == 0) {
  print('No arguments provided')
  quit()
} else {
  for (i in 1:length(args)){
    arg <- unlist(strsplit(args[i], '=', fixed=T))
    print(arg[1])
    if(arg[1] == 'input_dir') input_dir <- arg[2]
    if(arg[1] == 'output_dir') output_dir <- arg[2]
    if(arg[1] == 'regions') regions <- arg[2]
  }
}

if (input_dir == '') {
  print('No input dir provided')
  quit()
} else if (regions == '') {
  print('No region files provided')
  quit()
}

#### library load
library(GenomicRanges)

#### nameList generation
files <- list.files(input_dir)
tmp <- t(as.data.frame(strsplit(files, '.', fixed=T)))
nameList <- as.vector(tmp[,2])
prefix <- unique(tmp[,1])
suffix <- paste(unique(tmp[,3]), unique(tmp[,4]), sep='.')

#### regions load
anno <- readRDS(regions)
kind <- strsplit(strsplit(regions, '.', fixed=T)[[1]][1], '/', fixed=T)[[1]]
kind <- kind[length(kind)]

#### Functions
counts <- function(x)
{
  x <- as.vector(x)
  conC <- round(x[2] * x[1] / 100)
  conT <- round(x[3] * x[1] / 100)
  total <- conC + conT
  result <- c(conC, total)
  return(result)
}

averageMet <- function(x, met.matrix)
{
  if (length(x) == 0){return(c(0, 0))}
  else
  {
    tmp <- met.matrix[x,]
    tmp <- as.matrix(tmp[,5:7])
    tmp <- t(apply(tmp, 1, counts))
    tmp <- apply(tmp, 2, sum)
    methyl <- tmp[1]/tmp[2]
    return(c(methyl, tmp[1], tmp[2], length(x)))
  }
}

##   DNA Methylation Levels on Annotation
metLevelAnno <- function(nameList, anno, kind, input_dir, output_dir, prefix, suffix){
  for (i in 1:length(nameList)){
    file <- paste0(input_dir, '/', prefix, '.', nameList[i], '.', suffix)
    met.table <- read.table(file, header=TRUE)
    CpG.anno <- GRanges(seqnames=Rle(met.table$chr), ranges=IRanges(start=met.table$base, end=met.table$base))

    tmp <- as.list(findOverlaps(anno, CpG.anno))
    tmp <- lapply(tmp, averageMet, met.table)
    tmp <- as.data.frame(t(as.data.frame(tmp)) )
    rownames(tmp) <- names(anno)
    colnames(tmp) <- c('met', 'C_read', 'T_read', 'num')
    tmp <- tmp[tmp$num > 0,]
    saveRDS(tmp, paste0(output_dir, '/', kind, '.', nameList[i], '.rds'))
    print(paste0(nameList[i], "is finished!"))
  }
}
metLevelAnno(nameList, anno, kind, input_dir, output_dir, prefix, suffix)
