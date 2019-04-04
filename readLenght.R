library(EDASeq)

genes <- read.csv("genes.txt",sep=",",  header=TRUE)
lenghts <- getGeneLengthAndGCContent(genes[[1]],'hsa')
