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
  outdir <- "/scratch/Shares/dowell/dbnascent/deseq2/human_celltype_control_baseline_esc_hela"
  # Sample list file
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

# Generate column data to insert into dds object

colTable <- sample_list[-c(1:2)]
rownames(colTable) <- sample_list$sample_id
for (i in 1:ncol(colTable)) {
  if (colnames(colTable)[i] == compaxis) {
    compfact <- as.vector(unique(colTable[,i]))
  }
}

# Generate pairwise comparisons for comparison axis

compvals <- matrix(nrow = sum(1:(length(compfact)-1)), ncol = 2)
iter = 1
for (i in 1:(length(compfact)-1)) {
  for (j in (i+1):length(compfact)) {
    compvals[iter,1] = compfact[i]
    compvals[iter,2] = compfact[j]
    iter = iter + 1
  }
}

# Filter to most transcribed isoform and generate DESeq count matrix

countData = countData[order(countData$GeneID),]
n_occur <- data.frame(table(countData$GeneID))
to_filt <- n_occur[n_occur$Freq > 1,colnames(n_occur)]
countData_filt <- countData[countData$GeneID %in% n_occur[n_occur$Freq == 1,"Var1"], colnames(countData)]

for (i in 1:nrow(to_filt)) {
  to_comp <- countData[countData$GeneID == to_filt$Var1[i],colnames(countData)]
  sumrows <- apply(to_comp[,4:length(to_comp)],1, sum, na.rm=TRUE)
  countData_filt <- bind_rows(countData_filt,to_comp[which.max(sumrows), colnames(to_comp)])
}
rm(n_occur,to_filt)

countTable <- countData_filt[-c(1:3)]
dds <- def_dds(countTable,colTable)

# Calculate TPM values
# This can be done independently of DESeq2 with the count file, as well
#RPK <- (countTable)/(countData_filt$Length/1000)
#per_mil_scaling_factor <- colSums(RPK)/1000000
#TPM <- as.data.frame(t(t(data.matrix(RPK))/per_mil_scaling_factor))
#TPM$transcript_id <- countData_filt$TranscriptID
#TPM$gene_id <- countData_filt$GeneID
#write.table(as.data.frame(TPM),file=paste(outdir,"/all_samples_TPM.txt", sep=""),sep="\t",quote=FALSE,row.names = FALSE)


# Do differential expression analysis and split up pairwise comparisons
# Write out results, MA plots

dds <- DESeq(dds)

for (i in 1:nrow(compvals)) {
  res <- results(dds, contrast = c(compaxis,compvals[i,1],compvals[i,2]))
  res_labeled <- data.frame(res, GeneID = countData_filt$GeneID, TranscriptID = countData_filt$TranscriptID)
  
  write.table(as.data.frame(res_labeled),file=paste(outdir,"/deseq2_results_comp_",compvals[i,1],"_",compvals[i,2],".txt", sep = ""),sep="\t",quote=FALSE,row.names = FALSE)
  
  png(paste(outdir,"/deseq2_maplot_comp_",compvals[i,1],"_",compvals[i,2],".png", sep = ""),type='cairo')
  plotMA(res, alpha=0.05, main="DESeq2")
  dev.off()
  graphics.off()
}
