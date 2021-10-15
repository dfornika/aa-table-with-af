#!/usr/bin/env python3

import argparse

import pysam

def gen_aa_bcsq(consequence, gene, protein):
    aa = '-'.join(map(str, [gene, protein]))
    return aa


def gen_aa_ann(consequence, gene, protein):
    aa = ""
    if protein != '':
        aa = '-'.join(map(str, [gene, protein]))
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
            protein = bcsq[5]
            output.append(consequence)
            output.append(gene)
            output.append(protein)
            aa = gen_aa_bcsq(consequence, gene, protein)
            output.append(aa)
        elif 'ANN' in rec.info:
            ann = rec.info['ANN'][0].split('|')
            # print(ann)
            consequence = ann[1]
            gene = ann[3]
            protein = ann[10][2:]
            output.append(consequence)
            output.append(gene)
            output.append(protein)
            aa = gen_aa_ann(consequence, gene, protein)
            output.append(aa)
        else:
            output.append('')
        
        print('\t'.join(map(str, output)))
    

if __name__ == '__main__':
    main()
