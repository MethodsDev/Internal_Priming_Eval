#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import argparse
import logging
import subprocess


def main():

    parser = argparse.ArgumentParser(description="examine 3' read internal priming",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--bam", type=str, required=True, help="bam file, coordinate-sorted")
    parser.add_argument("--genome_fa", type=str, required=True, help="ref genome fasta file")
    parser.add_argument("--min_read_counts", type=int, default=2, help="min number of reads with same 3' end")
    parser.add_argument("--up_flank_dist", type=int, default=30, help="bases upstream of end")
    parser.add_argument("--down_flank_dist", type=int, default=30, help="bases downstream of end")
    
    args = parser.parse_args()

    input_bam_file = args.bam
    genome_fa_file = args.genome_fa
    min_count = args.min_read_counts
    up_flank_dist = args.up_flank_dist
    down_flank_dist = args.down_flank_dist


    seq_extract_config = { "genome_fa" : genome_fa_file,
                           "min_count" : min_count,
                           "up_flank_dist" : up_flank_dist,
                           "down_flank_dist" : down_flank_dist }
    
    
    samreader = pysam.AlignmentFile(
        input_bam_file, "rb"
    )

    if  ( (not 'SO' in samreader.header.as_dict()['HD']) or
        samreader.header.as_dict()['HD']['SO'] != 'coordinate') :
        raise RuntimeError("Error, file: {} must be coordinate sorted".format(input_bam_file))



    prime3_end_counter = defaultdict(int)


    prev_chrom = None

    for read in samreader.fetch():
        chrom = samreader.get_reference_name(read.reference_id)
        ref_start = read.reference_start
        ref_end = read.reference_end
        #read_align_rend = read.get_blocks()[-1][1]
        strand = '+' if read.is_forward else '-'

        transcript_end = ref_start if strand == '-' else ref_end

        if prev_chrom is not None and prev_chrom != chrom:
            dump_chrom_end_counts(prime3_end_counter, seq_extract_config)
            prime3_end_counter.clear()

        prev_chrom = chrom
        token = f"{chrom}^{transcript_end}^{strand}"
        prime3_end_counter[token] += 1

    dump_chrom_end_counts(prime3_end_counter, seq_extract_config) # get last one


    sys.exit(0)


def dump_chrom_end_counts(prime3_end_counter, seq_extract_config):
    for token, count in prime3_end_counter.items():
        if count >= seq_extract_config['min_count']:
            chrom, coord, strand = token.split("^")
            seqregion = extract_flank_seq(chrom, coord, strand, seq_extract_config)
            polyA_signal = check_for_polyA_signal(seqregion)
            if polyA_signal:
                polyA_signal_lc = polyA_signal.lower()
                seqregion = seqregion.replace(polyA_signal, polyA_signal_lc)
            
            print(f"{chrom}\t{coord}\t{strand}\t{count}\t{seqregion}")

    return



def extract_flank_seq(chrom, coord, strand, seq_extract_config):

    coord = int(coord)
    
    if strand == '+':
        lend_coord = coord - seq_extract_config['up_flank_dist']
        rend_coord = coord + seq_extract_config['down_flank_dist']
    else:
        lend_coord = coord - seq_extract_config['down_flank_dist']
        rend_coord = coord + seq_extract_config['up_flank_dist']

    token = f"{chrom}:{lend_coord}-{rend_coord}"

    genome_fa = seq_extract_config['genome_fa']
    
    cmd = f"samtools faidx {genome_fa} {token}"

    seq_string = "".join(subprocess.check_output(cmd, shell=True).decode().split("\n")[1:])
    seq_string = seq_string.upper()

    if strand == '-':
        seq_string = revcomp(seq_string)
    
    return seq_string



def revcomp(seq):
    # from: https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement
    


def check_for_polyA_signal(seq_region):

    pattern='(A[AT]TAAA)'
    m =  re.search(pattern, seq_region)
    if m:
        motif_match = m.group(1)
        return motif_match
    else:
        return None
    

if __name__=='__main__':
    main()


    
