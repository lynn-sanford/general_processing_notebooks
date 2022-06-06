#!/usr/bin/env python
# coding: utf-8

# transcript_cluster_ends_v1.0.py
# This script takes in refseq annotation bedfile of genes and adds 
# transcripts more than 100 bases different at 3' or 5' ends to 
# new GTF file

# Import libraries

import math
import numpy as np
import pandas as pd

# Define refseq input, output basename

infile = "/Users/lysa8537/data/annotation/hg38_copy/hg38_refseq_genenames_included.bed"
basename = "hg38_refseq"

# Dump input file into variable

with open(infile) as f:
    region_input = []
    for line in f:
        region = tuple(line.strip().split("\t"))
        region_input.append(region)

cols = [
    "chr",
    "start",
    "stop",
    "transcriptID",
    "misc",
    "strand",
    "start2",
    "stop2",
    "misc2",
    "numex",
    "exlengths",
    "exstarts",
    "geneID",
]

regions = pd.DataFrame(region_input, columns=cols)

# Filter out alternate chromosomes

filt_regions = regions[~regions["chr"].str.contains("_")]
filt_regions["fivep"] = np.where(filt_regions["strand"] == '+', filt_regions["start"], filt_regions["stop"])

# Collapse overlaps

filt_regions["start_round"] = np.ceil(filt_regions["start"].astype(float)/100)
filt_regions["stop_round"] = np.ceil(filt_regions["stop"].astype(float)/100)

unique_tx = filt_regions.sort_values(
    by=["geneID", "chr", "fivep"]
)
unique_tx = unique_tx.drop_duplicates(subset=["chr", "start_round", "stop_round", "strand", "geneID"])

for i in range(50):
    for j in range(50):
        unique_tx["temp_start_round"] = np.ceil((filt_regions["start"].astype(float) + i)/100)
        unique_tx["temp_stop_round"] = np.ceil((filt_regions["stop"].astype(float) + j)/100)
        unique_tx = unique_tx.drop_duplicates(subset=["chr", "start_round", "temp_stop_round", "strand", "geneID"])
        unique_tx = unique_tx.drop_duplicates(subset=["chr", "temp_start_round", "stop_round", "strand", "geneID"])
        unique_tx = unique_tx.drop_duplicates(subset=["chr", "temp_start_round", "temp_stop_round", "strand", "geneID"])

# Make columns for writing GTF

unique_tx["ident"] = "hg38_refseq"
unique_tx["type"] = "gene_length"
unique_tx["miscout"] = "0"
unique_tx["period"] = "."
unique_tx["geneidout"] = 'gene_id "'
unique_tx["txidout"] = '"; transcript_id "'
unique_tx["final_quote"] = '"'
unique_tx["desc_out"] = unique_tx[['geneidout',
                                   'geneID',
                                   'txidout',
                                   'transcriptID',
                                   'final_quote'
                                  ]].T.agg("".join)

# Export data

outcols = ['chr',
    'ident',
    'type',
    'start',
    'stop',
    'miscout',
    'strand',
    'period',
    'desc_out'
]

indir = infile.split("/")
outdir = "/".join(indir[0:-1])
outfile = "".join([outdir, "/", basename, "_diff53prime.gtf"])

unique_tx.to_csv(
    outfile,
    sep='\t',
    columns=outcols,
    header=False,
    index=False,
)
