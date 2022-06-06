import os
import sys
import math
import statistics

import argparse


def bedtools_slop(tfit, chrom, sample_name, outdir, scan_window=1000000):
    """Wrapper for running bedtools slop
    Parameters
    ----------
    tfit : str (path)
        path to tfit results bedfile

    chrom : str (path)
        path to genome chromosome sizes file

    sample_name : str
        sample name for output file

    outdir : str (path)
        path to output directory

    scan_window : int (default 1000000)
        half-width of total window for gene scanning
        i.e. distance from mu to scan

    Returns
    -------
    Output from 'bedtools slop'
    
    tfit_bed_slop : str (path)
        path to output bedfile
    """
    # set a base name for the output files

    base_name = sample_name
    tfit_bed_slop = str(outdir) + "/" + str(base_name) + "_tfit_for_gene_scan.bed"

    os.system(
        "bedtools slop -i {} -g {} -b {} > {}/{}_tfit_for_gene_scan.bed".format(
            tfit, chrom, scan_window, outdir, base_name
        )
    )
    
    return tfit_bed_slop


def read_bedfile(bed):
    """Read input BED files and return a list of coordinates

    Parameters
    ----------
    bed : str (path)
        path to bed file

    Returns
    -------
    bed_regions : list of lists
        coordinates from BED file

    """

    bed_regions = []

    with open(bed) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            # get individual regions from bed file
            regions = line.split("\t")

            # bed file from tfit results
            if len(regions) < 4:
                raise IndexError(
                    "Tfit BED file should have at least 4 columns. Double check input."
                )
            else:
                bed_regions.append(regions)

    return bed_regions


def mu_calc(tfit_regions):
    """Calculate midpoints of tfit regions

    Parameter
    ---------
    tfit_regions : list of lists
        regions parsed from tfit bedfile

    Returns
    -------
    mu_regions : list of lists
        regions of 1 bp corresponding to mu of the original input

    """
    mu_regions = []

    # calculate midpoint of region
    for region in tfit_regions:
        mu_coord = math.ceil(statistics.mean([int(region[1]),int(region[2])]))
        mu_regions.append([region[0], mu_coord, (mu_coord + 1), region[3]])
        
    return mu_regions


def parse_gtf(gtf, outdir, tss_window=1000):
    """load in gene gtf file and output TSS to bedfile

    Parameters
    ----------
    gtf : str (path)
        path to gene annotation file (in gtf format)

    outdir : str (path)
        output directory for gene bedfiles

    tss_window : int (default 1000)
        width of window around TSS to report

    Returns
    -------
    Bed file containing gene TSS regions

    tss_bed : str (path)
        path to output bedfile
    """

    # Load in gtf file and get 5' coordinate of gene
    tss_coord = []

    with open(gtf) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # get individual regions from bed file
            regions = line.split("\t")

            # check if gtf file
            if len(regions) != 9:
                raise IndexError(
                    "GTF file should have 9 columns. Double check input."
                )
            # Take first coordinate if on + strand
            elif regions[6] == "+":
                meta = regions[8].split("\"")
                region_to_add = [regions[0], (int(regions[3]) - int(math.ceil(tss_window/2))),
                                 (int(regions[3]) + int(math.ceil(tss_window/2))), meta[1]]
                if region_to_add not in tss_coord:
                    tss_coord.append([regions[0], (int(regions[3]) - int(math.ceil(tss_window/2))),
                                      (int(regions[3]) + int(math.ceil(tss_window/2))), meta[1]])
            # Take second coordinate if on - strand
            elif regions[6] == "-":
                meta = regions[8].split("\"")
                region_to_add = [regions[0], (int(regions[4]) - int(math.ceil(tss_window/2))),
                                 (int(regions[4]) + int(math.ceil(tss_window/2))), meta[1]]
                if region_to_add not in tss_coord:
                    tss_coord.append([regions[0], (int(regions[4]) - int(math.ceil(tss_window/2))),
                                      (int(regions[4]) + int(math.ceil(tss_window/2))), meta[1]])
            else:
                raise IndexError(
                    "GTF file should have strand info in 7th column. Double check input."
                )

    # Write out base 5' coordinate bedfile in bed4 format
    tss_bed = str(outdir) + "/gene_tss_" + str(tss_window) + "bp_window.bed"
    write_bedfile(tss_bed, tss_coord, 4)

    return tss_bed


def intersect_genes(tfit_bed_slop, tss_bed, sample_name, outdir, scan_window=1000000):
    """Wrapper for running bedtools intersect to intersect genes and TSS

    Parameters
    ----------
    tfit_bed_slop : str (path)
        path to bedfile of slop around mu for tfit regions

    tss_bed : str (path)
        path to bed file containing TSS regions

    sample_name : str
        basename for output files

    outdir : str (path)
        path to output directory

    scan_window : int (default 1000000)
        window used for slop file

    Returns
    -------
    Bedfile of intersections between windows around tfit regions
    and TSS of genes

    bed_intersect : str (path)
        path to intersect file
    """
    # set a base name for the output files and set paths
    base_name = str(sample_name)
    bed_intersect = (str(outdir) + "/" + base_name + 
                     "_tfit_gene_overlap_" + str(scan_window) + ".bed")

    os.system(
        "bedtools intersect -wa -wb -a {} -b {} > {}".format(
            tfit_bed_slop, tss_bed, bed_intersect
        )
    )

    return bed_intersect

                     
def unique_genes(intersect_regions):
    """Filter out any repeat gene overlaps from intersection. 

    Parameters
    ----------
    intersect_regions : list of lists
        regions derived from intersect bedfile

    Returns
    -------
    intersect_regions_filt : list of lists
        regions from intersect bedfile with gene repeats removed
    """
    intersect_regions_filt = []
    genes = []
    for region in intersect_regions:
        if region[7] in genes:
            continue
        else:
            intersect_regions_filt.append(region)
            genes.append(region[7])

    return intersect_regions_filt


def deseq_add(deseq_results, intersect_regions_filt, 
              significance=True, thresh=0.05):
    """Add DESeq2 results information to intersect region info.
    
    Parameters
    ----------
    deseq_results : str (path)
        path to deseq2 results file

    intersect_regions_filt : list of lists
        containing intersect regions with no gene repeats

    significance: boolean (default = True)
        whether to filter by significant padj

    thresh : float (default = 0.05)
        threshold for significance filtering

    Returns
    -------
    intersect_regions_filt with comparison, gene FC and padj appended
    """
    deseq = []
    genes = []
    file_split = deseq_results.split(".")[-2].split("_")
    comp = str(file_split[-2] + "_" + file_split[-1])

    with open(deseq_results) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # get individual gene entries
            gene_entries = line.split("\t")
            deseq.append(gene_entries)
            genes.append(gene_entries[6])

    for region in intersect_regions_filt:
        if region[7] not in genes:
            continue
        else:
            deseq_ext = deseq[genes.index(region[7])]
            if significance and (deseq_ext[5] == "NA"):
                continue
            elif significance and (float(deseq_ext[5]) >= thresh):
                continue
            else:
                region.append(comp)
                region.append(deseq_ext[7])
                region.append(deseq_ext[1])
                region.append(deseq_ext[5])

    return intersect_regions_filt


def write_bedfile(out_bedfile, region_list, columns=4):
    """Write output from the input list as a bedfile.
    Parameters
    ----------
    out_bedfile : str (file name)
        name of file where bed regions will be written

    region_list : list of lists
        containing regions in BED4 or BED5 format

    columns: int (default = 4)

    Returns
    -------
    writes 'out_bedfile' in BED4 or BED5 format
    """
    with open(out_bedfile, "w") as output:
        for region in region_list:
            if int(columns) == 4:
                output.write(
                    str(
                        str(region[0])
                        + "\t"
                        + str(region[1])
                        + "\t"
                        + str(region[2])
                        + "\t"
                        + str(region[3])
                        + "\n"
                    )
                )
            elif int(columns) == 5:
                output.write(
                    str(
                        str(region[0])
                        + "\t"
                        + str(region[1])
                        + "\t"
                        + str(region[2])
                        + "\t"
                        + str(region[3])
                        + "\t"
                        + str(region[4])
                        + "\n"
                    )
                )
            elif int(columns) == 12:
                output.write(
                    str(
                        str(region[0])
                        + "\t"
                        + str(region[1])
                        + "\t"
                        + str(region[2])
                        + "\t"
                        + str(region[3])
                        + "\t"
                        + str(region[4])
                        + "\t"
                        + str(region[5])
                        + "\t"
                        + str(region[6])
                        + "\t"
                        + str(region[7])
                        + "\t"
                        + str(region[8])
                        + "\t"
                        + str(region[9])
                        + "\t"
                        + str(region[10])
                        + "\t"
                        + str(region[11])
                        + "\n"
                    )
                )
            else:
                raise IndexError(
                    "Number of columns not supported for writing. Check input bedfile."
                )

        output.close()


def main(
    tfit, chrom, sample_name, output_dir, gtf, deseq_results, 
    tss_window=1000, scan_window=1000000, significance=True, thresh=0.05,
):
    """Process input prelim files"""

    # 1 : set a base name for the output files
    if sample_name:
        base_name = sample_name
    else:
        base = os.path.basename(tfit)
        base_name = base[: base.index(".")]

    # 2 : parse TSS regions
    tss_bed = parse_gtf(gtf, output_dir, tss_window)

    # 3 : load the tfit results bed file
    raw_tfit_results = read_bedfile(tfit)

    # 4: calculate mu and create slop bedfile
    tfit_mu = mu_calc(raw_tfit_results)
    mu_filename = output_dir + "/" + base_name + "_tfit_mu.bed"
    write_bedfile(mu_filename, tfit_mu, 4)
    comp_tfit_results = bedtools_slop(mu_filename, chrom, base_name, 
                                      output_dir, scan_window)

    # 5 : do intersection
    intersect_file = intersect_genes(comp_tfit_results, tss_bed, base_name, 
                                     output_dir, scan_window)
    
    # 6 : read in intersect file and filter for unique genes
    raw_intersect = read_bedfile(intersect_file)
    filtered_intersect = unique_genes(raw_intersect)

    # 7 : add deseq2 info
    deseq_intersect = deseq_add(deseq_results, filtered_intersect, 
                                significance, thresh)

    # 8 : write out results
    if significance:
        outfile = base_name + "_gene_scan_" + str(scan_window) + "_pfilt.bed"
    else:
        outfile = base_name + "_gene_scan_" + str(scan_window) + ".bed"

    write_bedfile(
        "{}/{}".format(output_dir, outfile),
        deseq_intersect,
        12,
    )


parser = argparse.ArgumentParser(description="Scan for genes around tfit results")

parser.add_argument(
    "-i", "--input", dest="tfit", help="Input tfit results BED file", metavar="FILE",
)
parser.add_argument(
    "-c", "--chrom", dest="chrom", help="Genome chromosome size file", metavar="FILE",
)

parser.add_argument(
    "-s", "--samp_id", dest="sample_name", help="Base name for sample", metavar="STR",
)

parser.add_argument(
    "-o", "--outdirectory", dest="output_dir", help="Directory for output", metavar="DIR",
)

parser.add_argument(
    "-g", "--gtf_file", dest="gtf", help="Gene annotion GTF file", metavar="FILE",
)

parser.add_argument(
    "-d", "--deseq_res", dest="deseq_results", help="Path to DESeq2 rsults", metavar="INT",
)

parser.add_argument(
    "-t",
    "--tss_window",
    dest="tss_window",
    type=int,
    default=1000,
    help="Width of window around TSS",
    metavar="INT",
)

parser.add_argument(
    "-w",
    "--scan_window",
    dest="scan_window",
    type=int,
    default=1000000,
    help="Distance from mu to scan for genes",
    metavar="INT",
)

parser.add_argument(
    "-p",
    "--pval",
    dest="thresh",
    type=float,
    default=0.05,
    help="Threshold for determining DESeq2 significance",
    metavar="FLOAT",
)

parser.add_argument(
    "--sig",
    dest="significance",
    action="store_true",
    help="Filter for significantly different genes",
)
parser.add_argument(
    "--no-sig",
    dest="significance",
    action="store_false",
    help="Report all genes, regardless of difference",
)

parser.set_defaults(significance=True)

args = parser.parse_args()

main(
    args.tfit,
    args.chrom,
    args.sample_name,
    args.output_dir,
    args.gtf,
    args.deseq_results,
    args.tss_window,
    args.scan_window,
    args.significance,
    args.thresh,
)
