#!/usr/bin/env python
import os
import argparse
from collections import Counter, defaultdict
from pprint import pprint

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Custom loader
from load_genotypes import load_genotypes

from Mappings import genotypeToHosts, species2group, manualMap

accToGenotype = load_genotypes()


def arg_parse() -> object:
    parser = argparse.ArgumentParser(description='HIV sequence analysis')
    parser.add_argument('-p1', type=str, default=None,
                        help='Path of ascension file 1')
    parser.add_argument('-p2', type=str, default=None,
                        help='Path of ascension file 2')
    return parser.parse_args()


# Plots a dict d object
# filename: the file name to write
# xName: what are you plotting (for title)
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
    plt.title(f'Frequency of {xName}', fontsize=20)

    # Print image
    name = f'./figures/{filename}-distribution.png'
    plt.savefig(name, bbox_inches='tight', dpi=300)
    plt.clf()


## -------------------------------------- Main method start --------------------------------------- ##


def main(args):
    # Manual mappings from isolate to host
    df = pd.read_csv("hosts.csv", sep=",", encoding="UTF-8")
    isolateLookupFailures = []

    def isolateHostLookup(isolate, acc):
        rowEntry = df.loc[df['GENBANK ACCESSION NUMBER'] == isolate]
        if len(rowEntry) != 0:
            return rowEntry['HOST'].iloc[0], True
        else:
            isolateLookupFailures.append((isolate, acc))
            return "NoIsolateLookupExists", False

    """ -------------- Parse list of records ----------------"""
    Entrez.email = 'andrewclchan211@vt.edu'
    Entrez.apikey = "62121718eb8c662846e0fbddb67a34480408"

    # If no accession path was passed in
    if not args.p1 or not args.p2:
        print("Error: no accession list found!")
        return
    else:
        print(f"Reading from file: {args.p1} as dataset 1")
        with open(args.p1, "r") as file:
            recordsIds1 = set(file.read().splitlines())

        print(f"Reading from file: {args.p2} as dataset 2")
        with open(args.p2, "r") as file:
            recordsIds2 = set(file.read().splitlines())

        in1Not2 = sorted(list(recordsIds1 - recordsIds2))
        in2Not1 = sorted(list(recordsIds2 - recordsIds1))

        print(f"In database 1, but not in 2: {in1Not2} of size {len(in1Not2)}")
        print("=================================================================")
        print(f"In database 2, but not in 1: {in2Not1} of size {len(in2Not1)}")




## -------------------------------------- Main method end --------------------------------------- ##

if __name__ == "__main__":
    args = arg_parse()
    print("Starting program")
    main(args)
