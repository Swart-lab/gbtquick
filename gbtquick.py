#!/usr/bin/env python3

import sys
import re
import argparse
import logging
import json
from collections import defaultdict
from statistics import mean,stdev,median

parser = argparse.ArgumentParser(description="""
    Report covstats and plot blobplots from (meta)genome assembly files without
    read-mapping, by parsing coverage info directly from headers
    """)
parser.add_argument("-a","--assembler", type=str, 
        help="Assembler used, valid options 'spades','flye','mapping'")
parser.add_argument("-f", "--fasta", type=str, 
        help="Assembly Fasta file")
parser.add_argument("-i", "--info", type=str, 
        help="assembly_info file from Flye")
parser.add_argument("--bam", type=str, 
        help="Mapping BAM file to calculate covstats")
parser.add_argument("-o", "--out", type=str, default="test", 
        help="Output filename prefix")
parser.add_argument("--dump", action="store_true",
        help="Dump data as json for troubleshooting")
args = parser.parse_args()

def fasta_to_gccount(filename):
    """Parse Fasta file and get GC base count per contig

    Parameters
    ----------
    filename : str
        Path to Fasta file

    Returns
    -------
    dict
        dict of ints (GC base counts) per contig, keyed by str (contig name)
    """
    # counter = 0
    gccount = defaultdict(int)
    with open(filename, "r") as fh:
        rname = ""
        for line in fh:
            line = line.rstrip() # Remove trailing whitespace
            if (re.match(">",line)):
                # counter += 1
                # if (counter % 50000 == 0):
                #     fhlog.write(str(counter)+ " entries parsed\n")
                rheader = line[1:] # Remove initial >
                # Strip anything after first whitespace
                rname = re.match("\S+", rheader)[0] 
            else:
                gccount[rname] += sum([line.count(char) for char in ["g","c","G","C"]])
                # TODO ambiguous IUPAC code for G/C
    return(gccount)


def parse_spades_assembly(assem):
    """Parse SPAdes assembly, get coverage from header and report covstats

    Parameters
    ----------
    assem : str
        Path to SPAdes scaffolds.fasta file
    
    Returns
    -------
    dict
        dict of dicts, keyed by str (statistic type), each dict contains stats
        values keyed by str (contig name)
    """
    covstats = defaultdict(dict)
    logging.info(f"Parsing SPAdes scaffolds assembly file {assem}")
    with open(assem) as fh:
        for line in fh:
            if re.match(r"^>", line):
                line = line.rstrip()
                rname = line[1:] # Strip > character from head
                rcap = re.match(r"NODE_(\d+)_length_(\d+)_cov_([\d.])", rname)
                if rcap:
                    covstats[rname]["Length"] = rcap.group(2)
                    covstats[rname]["Avg_fold"] = rcap.group(3)
                else:
                    logging.warn(f"Invalid SPAdes header format {rname}")
    logging.info(f"Parsing SPAdes scaffolds assembly file {assem} for GC content")
    gccount = fasta_to_gccount(assem)
    # Divide raw GC count by contig length to get GC frac
    for rname in covstats:
        if gccount[rname]:
            covstats[rname]["Ref_GC"] = float(gccount[rname] / covstats[rname]["Length"])
        else:
            covstats[rname]["Ref_GC"] = 0.0
    return(covstats)


def parse_flye_assembly(info, assem):
    """Parse Flye assembly, get coverage from header and report covstats

    Parameters
    ----------
    info : str
        Path to Flye assembly_info file
    assem : str
        Path to Flye assembly Fasta file

    Returns
    -------
    dict
        dict of dicts, keyed by str (statistic type), each dict contains stats
        values keyed by str (contig name)
    """
    covstats = defaultdict(dict)
    logging.info(f"Parsing Flye assembly_info file {info}")
    with open(info) as fh:
        for line in fh:
            if not re.match(r"^#", line): # skip comments/headers
                splitline = line.split(sep="\t")
                # seqname length cov circ repeat mult alt_group graph_path
                rname = splitline[0]
                covstats[rname]["Length"] = int(splitline[1])
                covstats[rname]["Avg_fold"] = float(splitline[2])
                covstats[rname]["Circular"] = str(splitline[3])
                covstats[rname]["Repeat"] = str(splitline[4])
    logging.info(f"Parsing Flye assembly Fasta file {assem} for GC content")
    gccount = fasta_to_gccount(assem)
    # Divide raw GC count by contig length to get GC frac
    for rname in covstats:
        if gccount[rname]:
            covstats[rname]["Ref_GC"] = float(gccount[rname] / covstats[rname]["Length"])
        else:
            covstats[rname]["Ref_GC"] = 0.0
    return(covstats)


def covstats_to_tsv(d, filename):
    """Convert covstats dict to TSV format

    Headers for all columns will be written
    """
    # scaffold Length Avg_fold Ref_GC Circular Repeat Cds_dens
    headings = ["Length", "Avg_fold", "Ref_GC", "Circular", "Repeat", "CDS_dens"]
    out = []
    heading_counter = defaultdict(int)
    for scaffold in d:
        outl = [scaffold]
        for heading in headings:
            if heading in d[scaffold]:
                heading_counter[heading] += 1
                outl.append(str(d[scaffold][heading]))
        out.append(outl)
    # Sanity check
    for heading in heading_counter:
        if heading_counter[heading] != len(d):
            logging.warn("Data lengths do not match number of scaffolds")
    # Write header
    header = ["Scaffold"]
    for heading in headings:
        if heading in heading_counter:
            header.append(heading)
    # Write file
    with open(filename, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for line in out:
            fh.write("\t".join(line) + "\n")


# main

if args.assembler == "flye":
    if not args.info:
        logging.warn("Flye assembly_info file not specified")
    covstats = parse_flye_assembly(args.info, args.fasta)
elif args.assembler == "spades":
    logging.log(f"Parsing SPAdes assembly file {args.fasta}")
    covstats = parse_spades_assembly(args.fasta)
else:
    logging.warn(f"Invalid assembler {args.assembler} specified")

if args.dump:
    logging.info("Dumping data to {args.out}.dump.json for troubleshooting")
    with open(f"{args.out}.dump.json","w") as fh:
        json.dump(covstats, fh, indent=4)
if covstats:
    covstats_to_tsv(covstats,f"{args.out}.covstats.tsv")
