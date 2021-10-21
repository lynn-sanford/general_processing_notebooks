# general_processing_notebooks
Mostly unrelated Jupyter notebooks for processing different aspects of seq and genome data

Each individual script has a summary of its function, but following is a basic description:

DESeq2_processing_for_intersect:
Takes DESeq2 results from non-gene regions and makes bedfiles for easy visualization of differential regions

TSS_bed_generate:
Creates bedfiles with windows around start codons input from refseq annotation files

add_refseq_genenames:
Adds gene names to refseq files that only have transcript numbers

make_genome_bedfile:
Converts a refseq gtf file to a bedfile

midpt_calc:
Within specific regions, calculates the midpt between the maxima from each half of the region

tpm_subset:
Subsets a large TPM or countfile based on an identifier list (currently written for subsetting TFs from all genes)

window_calc:
Outputs a bedfile with a specific window around region midpts of an input bedfile
