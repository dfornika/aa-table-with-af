#!/usr/bin/env python3

import argparse
import re

import pysam


def gen_protein_bcsq(consequence, protein):
    if consequence == 'synonymous':
        ref_aa = protein[-1]
        aa_change = ref_aa + protein
    elif consequence == 'missense':
        [ref_pos_aa, alt_pos_aa] = protein.split('>')
        aa_change = ref_pos_aa[-1] + alt_pos_aa
    elif re.search('deletion', consequence):
        aa_change = protein
    else:
        aa_change = protein

    return aa_change


def gen_aa(gene, aa_change):
    if aa_change != '':
        aa = '-'.join(map(str, [gene, aa_change]))
    else:
        aa = 'NA'
    return aa


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf')
    parser.add_argument('-s', '--sample_id')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    output_fields = [
        'sample',
        'chr',
        'pos',
        'ref',
        'alt',
        'depth',
        'alt_frequency',
        'consequence',
        'gene',
        'protein',
        'aa',
    ]

    print('\t'.join(output_fields))

    for rec in vcf.fetch():
        output = []
        if args.sample_id:
            sample_id = args.sample_id
        else:
            sample_id = rec.samples[0].name
        output.append(sample_id)
        output.append(rec.chrom)
        output.append(rec.pos)
        output.append(rec.ref)
        output.append(rec.alts[0])
        for s in rec.samples.items():
            output.append(s[1]['DP'])
        output.append(rec.info['AF'][0])
        if 'BCSQ' in rec.info:
            bcsq = rec.info['BCSQ'][0].split('|')
            consequence = bcsq[0]
            gene = bcsq[1]
            protein = gen_protein_bcsq(consequence, bcsq[5])
            output.append(consequence)
            output.append(gene)
            output.append(protein)
            aa = gen_aa(gene, protein)
            output.append(aa)
        elif 'ANN' in rec.info:
            ann = rec.info['ANN'][0].split('|')
            consequence = ann[1]
            gene = ann[3]
            protein = ann[10][2:]
            output.append(consequence)
            output.append(gene)
            output.append(protein)
            aa = gen_aa(gene, protein)
            output.append(aa)
        else:
            output.append('')
            output.append('')
            output.append('')
            output.append('NA')
        
        print('\t'.join(map(str, output)))
    

if __name__ == '__main__':
    main()
