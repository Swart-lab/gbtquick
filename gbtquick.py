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
    gccount = {}
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
    covstats = {}
    logging.info(f"Parsing Flye assembly_info file {info}")
    with open(info) as fh:
        for line in fh:
            if not re.match(r"^#", line): # skip comments/headers
                splitline = line.split(sep="\t")
                # seqname length cov circ repeat mult alt_group graph_path
                rname = splitline[0]
                covstats["Length"][rname] = float(splitline[1])
                covstats["Avg_fold"][rname] = float(splitline[2])
                covstats["Circular"][rname] = str(splitline[3])
                covstats["Repeat"][rname] = str(splitline[4])
    logging.info(f"Parsing Flye assembly Fasta file {assem} for GC content")
    gccount = fasta_to_gccount(filename)
    # Divide raw GC count by contig length to get GC frac
    for rname in covstats["Length"]:
        if gccount[rname]:
            covstats["Ref_GC"][rname] = float(gccount[rname] / covstats["Length"][rname])
        else:
            covstats["Ref_GC"][rname] = 0.0
    return(covstats)

# main

if args.assembler == "flye":
    # do flye things
    if not args.info:
        logging.warn("Flye assembly_info file not specified")
    covstats = parse_flye_assembly(args.info, args.fasta)

elif args.assembler == "spades":
    # do spades things
    logging.log(f"Parsing SPAdes assembly file {args.fasta}")
else:
    logging.warn(f"Invalid assembler {args.assembler} specified")
