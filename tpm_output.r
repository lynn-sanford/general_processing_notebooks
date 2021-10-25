library("DESeq2")
library("dplyr")

# User defined variables 
  # {"stranded","unstranded"}
  strand <- "stranded"
  # {"","5ptrunc_"}
  trunc <- ""
  # comparison axis
  compaxis <- "cell_type"
  # Output dir
#  outdir <- "C:/Users/Zarko/Desktop"
  outdir <- "/scratch/Shares/dowell/dbnascent/deseq2/human_celltype_control_baseline_esc_hela"
  # Sample list file
#  samplist <- "C:/Users/Zarko/Desktop/db_query_all_human.csv"
  samplist <- "/scratch/Shares/dowell/dbnascent/deseq2/human_celltype_control_baseline_esc_hela/db_query_human_control_baseline_ESC_HeLa.csv"
  # Function for defining DESeq2 design (likely need to edit design)
  def_dds <- function(countTable,colTable) {
    DESeqDataSetFromMatrix(countData = countTable, colData = colTable, design = ~ cell_type)
  }

# Take in list of samples and import all count data into one dataframe
  sample_list <- read.table(samplist,header=TRUE,sep =',')
  for (i in 1:nrow(sample_list)) {
#    new_ct <- read.table(paste("C:/Users/Zarko/Desktop/featurecounts_genes/",sample_list$sample_name[i],".sorted.",strand,".",trunc,"gene_counts.txt", sep=""),header=TRUE,sep='\t')
    new_ct <- read.table(paste("/Shares/dbnascent/",sample_list$paper_id[i],"/featurecounts_genes/",sample_list$sample_name[i],".sorted.",strand,".",trunc,"gene_counts.txt", sep=""),header=TRUE,sep='\t')
    sorted_ct <- new_ct[order(new_ct$TranscriptID),]
    if (i == 1) {
      countData <- sorted_ct
      colnames(countData) <- c(colnames(countData)[1:3], strsplit((colnames(countData)[4]),"\\.")[[1]][1])
      countData[,ncol(countData)] <- countData[,ncol(countData)] + 1
    } else {
      cname <- strsplit((colnames(sorted_ct)[4]),"\\.")[[1]][1]
      countData$newct <- sorted_ct[,ncol(sorted_ct)] + 1
      colnames(countData)[ncol(countData)] <- cname
    }
  }
rm(new_ct,sorted_ct)

# Filter to most transcribed isoform and generate DESeq count matrix

countTable <- countData[-c(1:3)]

# Calculate TPM values
# This can be done independently of DESeq2 with the count file, as well
RPK <- (countTable)/(countData$Length/1000)
per_mil_scaling_factor <- colSums(RPK)/1000000
TPM <- as.data.frame(t(t(data.matrix(RPK))/per_mil_scaling_factor))
TPM$transcript_id <- countData$TranscriptID
TPM$gene_id <- countData$GeneID
write.table(as.data.frame(TPM),file=paste(outdir,"/all_samples_TPM.txt", sep=""),sep="\t",quote=FALSE,row.names = FALSE)
write.table(as.data.frame(countData),file=paste(outdir,"/all_samples_countdata.txt", sep=""),sep="\t",quote=FALSE,row.names = FALSE)
