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
    parser.add_argument('-p', type=str, default=None,
                        help='Path of ascension files')
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
    maxUnitsInt = int(25 * (10 ** 3))

    # If no accession path was passed in
    if not args.p:
        print("Error: no accession list found!")
        return
    else:
        print(f"Reading from file: {args.p}")
        with open(args.p, "r") as file:
            recordsIds = file.read().splitlines()
            recordsIds = set(recordsIds[0:min(len(recordsIds), maxUnitsInt)])
    print(f"List of accession codes: {recordsIds}")

    ''' --------------- Query list of records for metadata --------------'''
    hostTypeMap = {}  # asccession : type
    notesWithoutHost = []
    proteinToNucAcc = {}

    print("Querying entries by ascension")
    idsRead = []
    ORF1 = []
    ORF2 = []
    ORF3 = []
    ORF4 = []
    unk = []

    tmpOut = [open(f'output/ORF{i}.fa', 'w') for i in range(1, 5)]

    # Get gb formatted data
    with Entrez.efetch(db="nuccore", id=",".join(recordsIds), rettype="gb", retmode="text") as handle:
        for seq_record in SeqIO.parse(handle, "gb"):  # Parse as genbank (gb) format
            idsRead.append(seq_record.annotations['accessions'][0])
            # Each entry has a sublist of sections, we only want the "source" which holds overall metadata
            for i, seq_feat in enumerate(seq_record.features):
                # Extraction of host features -------------------------------------
                if i == 0:
                    found = False
                    host = "Unknown"

                    if idsRead[-1] in manualMap:
                        host = manualMap[idsRead[-1]]
                        found = True

                    if 'host' in seq_feat.qualifiers:
                        host = seq_feat.qualifiers['host'][0]
                        found = True
                    if not found and 'lab_host' in seq_feat.qualifiers:
                        host = seq_feat.qualifiers['lab_host'][0]
                        found = True
                    if not found and 'isolate' in seq_feat.qualifiers:
                        key, status = isolateHostLookup(seq_feat.qualifiers['isolate'][0], idsRead[-1])
                        if status:  # If status code was true, we found it
                            host = key
                            found = True

                    if not found and 'isolation_source' in seq_feat.qualifiers:
                        key = seq_feat.qualifiers['isolation_source'][0]
                        if key in species2group:
                            host = key
                            found = True
                        else:
                            print(f"DEBUG: not found in species2group {key}")
                            isolateLookupFailures.append(key)

                    if not found and 'organism' in seq_feat.qualifiers:
                        for k in ("Avian", "Swine"):
                            if k.lower() in seq_feat.qualifiers['organism'][0].lower():
                                host = k
                                found = True
                                break
                    if idsRead[-1] in hostTypeMap:
                        print(f"WARN: duplicate sequences: {idsRead[-1]}")

                    # Assign 
                    if found:
                        if host in species2group and species2group[host] != "Unknown":
                            hostTypeMap[idsRead[-1]] = species2group[host]
                        else:
                            hostTypeMap[idsRead[-1]] = f"{host} (Unknown map)"
                            print("Unable to map to host:", idsRead[-1], seq_feat.qualifiers)
                            notesWithoutHost.append(seq_feat.qualifiers)
                    else:
                        hostTypeMap[idsRead[-1]] = host
                        print("Unable to map to host:", idsRead[-1], seq_feat.qualifiers)
                        notesWithoutHost.append(seq_feat.qualifiers)

                    # Only first feature, so always break

                # Extract ORFs
                if seq_feat.type == "CDS":
                    found = 0

                    for category in ("gene", "note", "product"):
                        if category in seq_feat.qualifiers:
                            if "ORF" in seq_feat.qualifiers[category][0].upper():
                                if "1" in seq_feat.qualifiers[category][0]:
                                    ORF1.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 1
                                elif "2" in seq_feat.qualifiers[category][0]:
                                    ORF2.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 2
                                elif "3" in seq_feat.qualifiers[category][0]:
                                    ORF3.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 3
                                elif "4" in seq_feat.qualifiers[category][0]:
                                    ORF4.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 4
                                if found:
                                    break
                            elif "protein" in seq_feat.qualifiers[category][0].lower():
                                if "poly" in seq_feat.qualifiers[category][0].lower() or (
                                        "non" and "structural" in seq_feat.qualifiers[category][0].lower()):
                                    ORF1.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 1
                                elif "capsid" in seq_feat.qualifiers[category][0].lower():
                                    ORF2.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 2
                                elif "hypothetical" in seq_feat.qualifiers[category][0].lower():
                                    ORF3.append(seq_feat.qualifiers['protein_id'][0])
                                    found = 3
                                if found:
                                    break

                    if not found:
                        unk.append(seq_feat.qualifiers)
                    else:
                        # Formatting features
                        proteinToNucAcc[seq_feat.qualifiers['protein_id'][0]] = idsRead[-1]
                        translation = seq_feat.qualifiers['translation'][0]
                        genotype = "unknown" if idsRead[-1] not in accToGenotype else accToGenotype[idsRead[-1]]
                        feature_sequence = SeqRecord(Seq(translation, ), id=seq_feat.qualifiers['protein_id'][0],
                                                     description="|ORF{}|{}|{}".format(found, hostTypeMap[idsRead[-1]],
                                                                                       genotype), )
                        # Write to file
                        SeqIO.write(feature_sequence, tmpOut[found - 1], 'fasta')

    map(lambda x: x.close(), tmpOut)  # Close all

    print(len(ORF1))
    print(len(ORF2))
    print(len(ORF3))
    print(len(ORF4))
    print(len(unk))

    hostTypeCounter = Counter(x for x in hostTypeMap.values())

    # Verbose print outs and matplotlib:
    # How many was entrex able to query from genbank
    print("Unable to querry the following accession codes: ", set(recordsIds) - set(idsRead),
          " Successful api querry: {}/{}".format(len(idsRead), len(recordsIds)))
    print(f'\n Host counts:  {hostTypeCounter} with size:', sum(x for x in hostTypeCounter.values()))
    print(f'\n\n DEBUG: Isolate lookup failures: {isolateLookupFailures} with size:', len(isolateLookupFailures))
    print(f'\n\n DEBUG: Unknown entries: {notesWithoutHost} with size:', len(notesWithoutHost))
    plotCounter(hostTypeCounter,  (os.path.basename(args.p) if args.p else "None") + "-post-processing-host", "host (processed)")

    hostGenotypeMiss = Counter()
    hostGenotypeMatch = Counter()

    unmappedHosts = set(x for x in hostTypeMap.values()) - set(e for x in genotypeToHosts.values() for e in x)
    unmappedHits = Counter()
    perGenotypeHostDist = defaultdict(Counter)  # HEVi : Counter(host)

    # See how many genotypes match with mapping
    for acc, refHost in hostTypeMap.items():
        if refHost in unmappedHosts:
            unmappedHits[refHost] += 1

        if acc not in accToGenotype or accToGenotype[acc] not in genotypeToHosts or refHost in unmappedHosts:
            print("Warn: no scrapping/mapping/genotypeing for ", acc)
        else:
            genotype = accToGenotype[acc]
            if refHost not in genotypeToHosts[genotype]:
                hostGenotypeMiss[genotype] += 1
                print("miss", refHost)
            else:
                perGenotypeHostDist[genotype][refHost] += 1
                hostGenotypeMatch[genotype] += 1

    print(f"PerGenotypeHostDistribution: {perGenotypeHostDist}")
    for genotype, fMap in perGenotypeHostDist.items():
        plotCounter(fMap, f"{genotype}-Host", f"{genotype}")

    # write data to csv

    df = pd.DataFrame.from_dict(perGenotypeHostDist)
    df.to_csv(r'genotypeHostDist.csv', index=False, header=True)


    missTotal = sum(hostGenotypeMiss.values())
    hitTotal = sum(hostGenotypeMatch.values())
    print(f"host type miss: {hostGenotypeMiss} total: {missTotal}")
    print(f"host type match: {hostGenotypeMatch} total: {hitTotal}")
    print(
        f"Valid comparisons (valid scraped hosts and valid calculated genotypes): {missTotal + hitTotal} / {len(idsRead)}")


### -------------------------------------- Main method end --------------------------------------- ##

if __name__ == "__main__":
    args = arg_parse()
    print("Starting program")
    main(args)
