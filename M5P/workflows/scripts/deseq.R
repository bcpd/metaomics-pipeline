#!/usr/bin/env Rscript

# Levy-Booth, Dimitriu
# For use with snakemake

#Load Libraries
##########################################################################################################
write("Loading DESeq2", stdout())
suppressMessages(library("DESeq2"))


parallel <- FALSE

#Load data
##########################################################################################################

for(count_table in snakemake@input[["counts"]]) {
  counts <- read.table(count_table, header=TRUE, row.names = 1, check.names=FALSE)
  coldata <- read.table(snakemake@params[["data"]], header=TRUE, check.names=FALSE)
  #all_conditions <- snakemake@params[["contrasts"]]
  #all_conditions <- snakemake@params[["design"]]

#Format data
##########################################################################################################
  counts <- as.matrix(counts); mode(counts) <- "integer"
  counts[is.na(counts)] <- 0

  all_conditions <- unique(coldata$treatment)

  cont_out <- data.frame(t(combn(all_conditions, 2)))
  #cont_out <- strsplit(as.character(all_conditions),'_')
  #cont_out <- data.frame(do.call(rbind, cont_out))

#Init DESeq2
##########################################################################################################
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ treatment)

  # remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  # normalization and preprocessing
  write("Running DESeq2", stdout())
  dds <- DESeq(dds, parallel=parallel, quiet = TRUE)

  for(i in 1:nrow(cont_out)) {
    c1 <- as.character(cont_out$X1[i])
    c2 <- as.character(cont_out$X2[i])
    #write(" ", stdout())
    write(paste("DeSeq2: Comparing ", c1, " (+) and ", c2, " (-)", sep = ""), stdout())

    contrast <- c("Condition", c1, c2)
    res <- results(dds, contrast=contrast, parallel=parallel)
    res <- suppressMessages(lfcShrink(dds, contrast=contrast, res=res))
    res <- res[order(res$log2FoldChange),]

    write(paste("DeSeq2: writing output:",snakemake@output[["tables"]][i]), stdout())

    #build filename
    table_file <- basename(snakemake@input[["counts"]])
    table_file <-  gsub("counts", paste(c1, c2, sep = "_"), table_file)
    #write.table(as.data.frame(res), file=paste(snakemake@params[["output_directory"]], "/tables/", table_file, sep = ""))
    write.table(as.data.frame(res), file=snakemake@output[["tables"]][i])
  }
}
