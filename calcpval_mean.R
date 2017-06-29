pvalmean <- function(infilenames, outfilename) {
	curtable <- read.table(infilenames[1], header=T)
	ptable <- curtable[,1:2]
	for(i in 2:length(infilenames)){
	      cat(i, "\n")
	      curtable <- read.table(infilenames[i], header=T)
      	      ptable <- merge(ptable, curtable[,1:2], by="gene")
	}

	curtable$pval <- rowMeans(ptable[2:ncol(ptable)])
	curtable$padj <- p.adjust(curtable$pval, "BH")

	write.table(curtable, file=outfilename, row.names=F, sep="\t", quote=F)
}
