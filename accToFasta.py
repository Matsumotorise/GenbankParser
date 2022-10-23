#!/usr/bin/env python
import argparse

from Bio import SeqIO, Entrez

import matplotlib.pyplot as plt
import os




def arg_parse():
    parser = argparse.ArgumentParser(description='HIV sequence analysis')
    parser.add_argument('-p', type=str, default=None,
                        help='Path of ascension files')
    return parser.parse_args()


# d: Counter(), p_name:
def plotCounter(d, filename, xName="genotype"):
    with plt.style.context("seaborn"):
        fig = plt.figure(1, [16, 9])
        # plt.rcParams.update({'font.size': 25})
        tmp = sorted(zip(d.keys(), d.values()), key=lambda x: -x[1])
        keys = [x[0] for x in tmp]
        values = [y[1] for y in tmp]
        plt.bar(keys, values)
        plt.xticks(
            rotation=45,
            horizontalalignment='right',
            # fontweight='heavy',
            # fontsize='small'
        )

    plt.xlabel(xName, fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.title(f'Frequency of each {xName} based on genbank notes', fontsize=20)

    # Print image
    name = f'{filename}-distribution.png'
    plt.savefig(name, bbox_inches='tight', dpi=300)
    plt.clf()


### -------------------------------------- Main method start --------------------------------------- ##


def main(args):
    ''' -------------- Parse list of records ----------------'''
    Entrez.email = 'andrewclchan211@vt.edu'
    Entrez.apikey = ""

    maxUnitsInt = int(25 * (10 ** 3))

    outFile = open(os.path.splitext(args.p)[0] + '.fasta', "w")

    # If no accession path was passed in
    if not args.p:
        return
    else:
        print(f"Reading from file: {args.p}")
        with open(args.p, "r") as file:
            recordsIds = file.read().splitlines()
            recordsIds = recordsIds[0:min(len(recordsIds), maxUnitsInt)]
    print(recordsIds)
    recordsIds = set(recordsIds)

    with Entrez.efetch(db="Protein", id=",".join(list(recordsIds)), rettype="gb", retmode="text") as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            seq_record.id = seq_record.annotations['accessions'][0]
            if seq_record.id not in recordsIds:
                print(f"WARN invalid seq acc: {seq_record.name}")
                print(seq_record)
            seq_record.description = ""

            SeqIO.write(seq_record, outFile, 'fasta')

    outFile.close()


### -------------------------------------- Main method end --------------------------------------- ##

if __name__ == "__main__":
    args = arg_parse()
    print("Starting program")
    main(args)
