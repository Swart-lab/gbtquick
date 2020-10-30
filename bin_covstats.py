#!/usr/bin/env python3

import argparse
import re
import logging

logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser(description="""
    From covstats.tsv file produced by gbtquick, filter Fasta file by criteria
    of length, GC, CDS density.
    """)
# input files
parser.add_argument("-f", "--fasta", type=str, 
    help="Assembly Fasta file")
parser.add_argument("--covstats", type=str,
    help="Covstats file produced by gbtquick")
# parameters
parser.add_argument("--min_length", type=int, default=None,
    help="Minimum contig length (bp), exclusive")
parser.add_argument("--max_length", type=int, default=None,
    help="Maximum contig length (bp), exclusive")
parser.add_argument("--min_gc", type=float, default=0.0,
    help="Minimum gc content, between 0 and 1")
parser.add_argument("--max_gc", type=float, default=1.0,
    help="Maximum gc content, between 0 and 1")
parser.add_argument("--min_cds_dens", type=float, default=None,
    help="Minimum cds density, between 0 and 1")
parser.add_argument("--max_cds_dens", type=float, default=None,
    help="Maximum cds density, between 0 and 1")
# output params
parser.add_argument("-o", "--out", type=str, default="test", 
    help="Output filename prefix")
parser.add_argument("--out_fmt", type=str, default="list",
    help="""
    Format for output file, either 'list' or 'fasta'. list output can be used
    with faSomeRecords from kentutils, for example, which is faster than this
    python implementation.
    """)
args = parser.parse_args()

def parse_covstats_tsv(filename):
    out = []
    fh = open(filename, "r")
    # get fields from header row
    headers = fh.readline().rstrip().split("\t")
    # process subsequent lines
    for line in fh:
        linesplit = line.rstrip().split("\t")
        out.append({headers[i]:linesplit[i] for i in range(len(headers))})
    fh.close()
    return(out)

def filter_field(recs, field, minval, maxval):
    out = []
    for rec in recs:
        if field in rec:
            if float(rec[field]) > minval and float(rec[field]) < maxval:
                out.append(rec)
        else:
            print(f"Field {field} not found in file!")
    return(out)

def filter_fasta(infile, outfile, headers):
    switch = False
    fh = open(infile, "r")
    fh_out = open(outfile, "w")
    counter = 0
    for line in fh:
        headmatch = re.match(r">(\S+)", line) # stop at first whitespace, behavior of input
        if headmatch:
            counter += 1
            if counter%1000 == 0:
                logging.info(f"... {str(counter)} Fasta entries processed")
            # check if header in list, turn on switch
            # else turn off switch
            if headmatch.group(1) in headers:
                switch=True
                fh_out.write(line)
            else:
                switch=False
        else:
            if switch:
                fh_out.write(line)
    fh.close()
    fh_out.close()

if __name__ == "__main__":
    recs = parse_covstats_tsv(args.covstats)
    logging.info(f"Total of {len(recs)} entries in covstats file {args.covstats}")
    if args.min_length or args.max_length:
        min_length = args.min_length if args.min_length else 0
        max_length = args.max_length if args.max_length else float('inf')
        logging.info(f"Filter: Length min {str(min_length)} max {str(max_length)}")
        recs = filter_field(recs, "Length", min_length, max_length)
        logging.info(f"After filter: {len(recs)} entries")
    if args.min_gc or args.max_gc:
        min_gc = args.min_gc if args.min_gc else 0.0
        max_gc = args.max_gc if args.max_gc else 1.0
        logging.info(f"Filter: Ref_GC min {str(min_gc)} max {str(max_gc)}")
        recs = filter_field(recs, "Ref_GC", min_gc, max_gc)
        logging.info(f"After filter: {len(recs)} entries")
    if args.min_cds_dens or args.max_cds_dens:
        min_cds_dens = args.min_cds_dens if args.min_cds_dens else 0.0
        max_cds_dens = args.max_cds_dens if args.max_cds_dens else 0.0
        logging.info(f"Filter: CDS_dens min {str(min_cds_dens)} max {str(max_cds_dens)}")
        recs = filter_field(recs, "CDS_dens", min_cds_dens, max_cds_dens)
        logging.info(f"After filter: {len(recs)} entries")
    headers = [rec['Scaffold'] for rec in recs]
    if args.out_fmt == "list":
        # report list of headers
        logging.info(f"Writing list of headers to {args.out}.filter.list...")
        with open(args.out+".filter.list", "w") as fh:
            for rec in headers:
                fh.write(rec + "\n")
    elif args.out_fmt == "fasta":
        logging.info(f"Filtering Fasta file {args.fasta}...")
        logging.info(f"Writing output to {args.out}.filter.fasta...")
        filter_fasta(args.fasta, args.out+".filter.fasta", headers)

