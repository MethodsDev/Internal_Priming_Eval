#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import argparse
import logging
import subprocess
import gzip


A_FRAC_WIN_SIZE = 20
KNOWN_POLYA_SEARCH_DIST = 50
UP_FLANK_DIST_DEFAULT = 30
DOWN_FLANK_DIST_DEFAULT = 30


POLYA_MOTIF_pattern='(A[AT]TAAA)'

SQANTI3_POLYA_MOTIFS = [x.upper() for x in  [
    # aataaa  # captured by regex above
    # attaaa  # captured by regex above
    'agtaaa',
    'tataaa',
    'cataaa',
    'gataaa',
    'aatata',
    'aataca',
    'aataga',
    'aaaaag',
    'actaaa',
    'aagaaa',
    'aatgaa',
    'tttaaa',
    'aaaaca',
    'ggggct' ] ]
# polyA motifs likely from here: https://polyasite.unibas.ch/atlas


DATADIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../data")



def main():

    parser = argparse.ArgumentParser(description="examine 3' read internal priming",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--bam", type=str, required=True, help="bam file, coordinate-sorted")
    parser.add_argument("--genome_fa", type=str, required=True, help="ref genome fasta file")
    parser.add_argument("--min_read_counts", type=int, default=2, help="min number of reads with same 3' end")
    parser.add_argument("--up_flank_dist", type=int, default=UP_FLANK_DIST_DEFAULT, help="bases upstream of end")
    parser.add_argument("--down_flank_dist", type=int, default=DOWN_FLANK_DIST_DEFAULT, help="bases downstream of end")
    
    args = parser.parse_args()

    input_bam_file = args.bam
    genome_fa_file = args.genome_fa
    min_count = args.min_read_counts
    up_flank_dist = args.up_flank_dist
    down_flank_dist = args.down_flank_dist


    genome_seq_lengths = parse_genome_seq_lengths(f"{genome_fa_file}.fai")

    known_polyA_sites = parse_known_polyA_sites()
    

    seq_extract_config = { "genome_fa" : genome_fa_file,
                           "min_count" : min_count,
                           "up_flank_dist" : up_flank_dist,
                           "down_flank_dist" : down_flank_dist,
                           "genome_seq_lengths" : genome_seq_lengths,
                           "known_polyA_sites" : known_polyA_sites
    }
    
    
    
    samreader = pysam.AlignmentFile(
        input_bam_file, "rb"
    )

    if  ( (not 'SO' in samreader.header.as_dict()['HD']) or
        samreader.header.as_dict()['HD']['SO'] != 'coordinate') :
        raise RuntimeError("Error, file: {} must be coordinate sorted".format(input_bam_file))



    prime3_end_counter = defaultdict(int)

    # print header for tabular report
    print("\t".join(["Chromosome",
                     "Coordinate",
                     "Strand",
                     "Num_reads",
                     "SeqRegion",
                     "EndCountA",
                     "EndFracA",
                     "PolyA_motif",
                     "AtOrNearKnownPolyA"]))
    
    
    prev_chrom = None

    for read in samreader.fetch():
        chrom = samreader.get_reference_name(read.reference_id)

        if chrom not in genome_seq_lengths:
            continue # skip those not in our genome fasta file.
        
        ref_start = read.reference_start
        ref_end = read.reference_end
        #read_align_rend = read.get_blocks()[-1][1]
        strand = '+' if read.is_forward else '-'

        transcript_end = ref_start if strand == '-' else ref_end

        if prev_chrom is not None and prev_chrom != chrom:
            dump_chrom_end_counts(prev_chrom, prime3_end_counter, seq_extract_config)
            prime3_end_counter.clear()

        prev_chrom = chrom
        token = f"{chrom}^{transcript_end}^{strand}"
        prime3_end_counter[token] += 1

    dump_chrom_end_counts(prev_chrom, prime3_end_counter, seq_extract_config) # get last one


    sys.exit(0)


def dump_chrom_end_counts(chromosome, prime3_end_counter, seq_extract_config):

    polyA_sites = seq_extract_config['known_polyA_sites'][chromosome]
    
    for token, count in prime3_end_counter.items():
        if count >= seq_extract_config['min_count']:
            chrom, coord, strand = token.split("^")
            coord = int(coord)

            # exclude earlier out-of-range polyA sites
            while len(polyA_sites)>0 and coord - polyA_sites[0]['coord'] > KNOWN_POLYA_SEARCH_DIST:
                polyA_sites.pop(0)

            at_known_polyA_site = False
            for polyA_site in polyA_sites:
                if (abs(coord - polyA_site['coord']) <= KNOWN_POLYA_SEARCH_DIST and
                    strand == polyA_site['strand']):
                    at_known_polyA_site = True
                    break
                if polyA_site['coord'] - coord > KNOWN_POLYA_SEARCH_DIST:
                    break
                            
            seqregion, rel_end_coord = extract_flank_seq(chrom, coord, strand, seq_extract_config)

            count_A, frac_A = compute_frac_A(seqregion, rel_end_coord, A_FRAC_WIN_SIZE)

            polyA_signal = check_for_polyA_signal(seqregion)
            if polyA_signal:
                polyA_signal_lc = polyA_signal.lower()
                seqregion = seqregion.replace(polyA_signal, polyA_signal_lc)
            
            print(f"{chrom}\t{coord}\t{strand}\t{count}\t{seqregion}\t{count_A}\t{frac_A:.3f}\t{polyA_signal}\t{at_known_polyA_site}")

    return



def extract_flank_seq(chrom, coord, strand, seq_extract_config):

    coord = int(coord)
    
    if strand == '+':
        lend_coord = coord - seq_extract_config['up_flank_dist']
        rend_coord = coord + seq_extract_config['down_flank_dist']
    else:
        lend_coord = coord - seq_extract_config['down_flank_dist']
        rend_coord = coord + seq_extract_config['up_flank_dist']

    lend_coord = max(1, lend_coord)
    rend_coord = min(rend_coord, seq_extract_config['genome_seq_lengths'][chrom])

    rel_end_coord = coord - lend_coord
    
    token = f"{chrom}:{lend_coord}-{rend_coord}"

    genome_fa = seq_extract_config['genome_fa']
    
    cmd = f"samtools faidx {genome_fa} {token}"

    seq_string = "".join(subprocess.check_output(cmd, shell=True).decode().split("\n")[1:])
    seq_string = seq_string.upper()

    if strand == '-':
        seq_string = revcomp(seq_string)
        rel_end_coord = len(seq_string) - rel_end_coord + 1
        
    return seq_string, rel_end_coord



def revcomp(seq):
    # from: https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement
    


def check_for_polyA_signal(seq_region):

    #pattern='(A[AT]TAAA)'
    m =  re.search(POLYA_MOTIF_pattern, seq_region)
    if m:
        motif_match = m.group(1)
        return motif_match
    else:
        for motif in SQANTI3_POLYA_MOTIFS:
            if seq_region.find(motif) >= 0:
                return motif
        return None
    


def parse_genome_seq_lengths(genome_fai_file):

    genome_seq_lengths = dict()

    with open(genome_fai_file, "rt") as fh:
        for line in fh:
            vals = line.split("\t")
            chrom = vals[0]
            seqlen = vals[1]
            genome_seq_lengths[chrom] = int(seqlen)

    return genome_seq_lengths


def compute_frac_A(seqregion, end_coord, win_size):
    
    lend_coord = end_coord
    rend_coord = min(end_coord+win_size, len(seqregion))

    seqregion = seqregion[lend_coord-1 : rend_coord]

    if len(seqregion) == 0:
        return "NA", "NA"
    
    A_count = 0
    for char in seqregion:
        if char == "A":
            A_count += 1
    frac_A = A_count / len(seqregion)

    return A_count, frac_A


def parse_known_polyA_sites():

    polyA_sites_file = os.path.join(DATADIR, "atlas.clusters.2.0.GRCh38.96.bed.gz")

    known_polyA_sites = defaultdict(list)

    with gzip.open(polyA_sites_file, "rt") as fh:
        for line in fh:
            vals = line.split("\t")
            token = vals[3]
            chrom, coord, strand = token.split(":")
            chrom = f"chr{chrom}"
            known_polyA_sites[chrom].append( { 'coord' : int(coord),
                                               'strand' : strand } )


    # ensure sorted
    chroms = list(known_polyA_sites.keys())
    for chrom in chroms:
        known_polyA_sites[chrom] = sorted(known_polyA_sites[chrom], key=lambda x: x['coord'])

    return known_polyA_sites

        
if __name__=='__main__':
    main()


    
