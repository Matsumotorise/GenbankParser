#!/usr/bin/env python
import os
import regex as re
import sys
import json
import argparse
import collections

from Bio import SeqIO, Entrez

import pandas as pd

import matplotlib.pyplot as plt

from pprint import pprint

def arg_parse():
    parser = argparse.ArgumentParser(description='HIV sequence analysis')
    parser.add_argument('-p', type=str, default=None,
                        help='Path of ascension files')
    return parser.parse_args()



# d: Counter(), p_name: 
def plotCounter(d, filename, xName="genotype"):
    with plt.style.context("seaborn"):
        fig = plt.figure(1, [16, 9])
        #plt.rcParams.update({'font.size': 25})
        tmp = sorted(zip(d.keys(), d.values()), key=lambda x: -x[1])
        keys = [x[0] for x in tmp]
        values = [y[1] for y in tmp]
        plt.bar(keys, values)
        plt.xticks(
            rotation=45,
            horizontalalignment='right',
            #fontweight='heavy',
            #fontsize='small'
        )

    plt.xlabel(xName, fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.title(f'Frequency of each {xName} based on genbank notes',fontsize=20)

    # Print image
    name = f'{filename}-distribution.png'
    plt.savefig(name, bbox_inches='tight', dpi=300)
    plt.clf()



### -------------------------------------- Main method start --------------------------------------- ##


def main(args):
    ''' -------------- Parse list of records ----------------'''
    Entrez.email='andrewclchan211@vt.edu'
    Entrez.apikey= ""

    maxUnitsInt = int(25*(10**3))
    maxUnits = str(maxUnitsInt)


    df = pd.read_csv("hosts.csv", sep=",", encoding="UTF-8")

    isolateLookupFailures = []
    def isolateHostLookup(isolate):
        rowEntry = df.loc[df['GENBANK ACCESSION NUMBER'] == isolate]
        if len(rowEntry) != 0:
            return rowEntry['HOST'].iloc[0]
        else:
            isolateLookupFailures.append(isolate)
            return "NoIsolateLookupExists"



    # If no accention path was passed in 
    if not args.p:
        virus = "Orthohepevirus A[Organism]"
        filterQuery = "AND (complete[Properties] OR complete[Title])"
        # Read records with complete only and assign "acc" accention codes to records
        with Entrez.esearch(db="protein", term=f"{virus} {filterQuery}", idtype="acc", retmax=maxUnits) as handle:
            records = Entrez.read(handle)

        assert records
        assert records['IdList']
        # Print out head of record list
        print(records['IdList'][0:min(maxUnitsInt, 10)], len(records['IdList']))
        # hashmap recordIds
        recordsIds = set(records['IdList'])
    else:
        print(f"Reading from file: {args.p}")
        with open(args.p, "r") as file:
            recordsIds = file.read().splitlines()
            recordsIds = recordsIds[0:min(len(recordsIds), maxUnitsInt)]
    print(recordsIds)




    ''' --------------- Query list of records for metadata --------------'''
    regex = r"(?:g|G)enotype(?:\s|:|-)*(\(.*\)|.*?)(?:$|\s|;|/)" # Capture group to search for in notes if G/genotype exists
    genoTypeDict = collections.Counter() # Make frequency dictionaries
    hostTypeDict = collections.Counter()
    notesWithoutGenotype = [] # For debuggging in case our search is missing something
    notesWithoutHost = []


    print("Querying entries by ascension")
    idsRead = []
    # Get gb formatted data
    with Entrez.efetch(db="Protein", id=",".join(recordsIds), rettype="gb", retmode="text") as handle:
        for seq_record in SeqIO.parse(handle, "gb"): # Parse as genbank (gb) format
            idsRead.append(seq_record.id)
            for i, seq_feat in enumerate(seq_record.features): # Each entry has a sublist of sections, we only want the "source" which holds overall metadata
                # Extraction of genotype if it exists ------------------------------------------
                if 'serotype' in seq_feat.qualifiers:
                    print(f'DEBUG serotype field found: {seq_feat.qualifiers}')

                # If there's a note, and if the note contains 'genotype', try to extract genotype
                if 'note' in seq_feat.qualifiers:
                    matches = re.findall(regex, seq_feat.qualifiers['note'][0], overlapped=True)
                    if len(matches) >= 1:
                        genoTypeDict[matches[-1]] += 1
                        if any('genotype' in m or '' == m for m in matches):
                            print(f'DEBUG:', seq_feat.qualifiers['note'][0], matches)
                        if len(matches) > 1:
                            print("warn: more than one regex match found: " + str(matches) +  seq_feat.qualifiers['note'][0])
                    else:
                        notesWithoutGenotype.append(seq_feat.qualifiers['note'][0])
                        genoTypeDict['None'] += 1
                elif 'genotype' in seq_feat.qualifiers:
                    print(f'DEBUG genotype field found: {seq_feat.qualifiers}')
                else:
                    genoTypeDict['None'] += 1

                # Extraction of host featuers -------------------------------------
                if 'host' in seq_feat.qualifiers:
                    hostTypeDict[seq_feat.qualifiers['host'][0]] += 1
                elif 'lab_host' in seq_feat.qualifiers:
                    hostTypeDict[seq_feat.qualifiers['lab_host'][0]] += 1
                elif 'isolate' in seq_feat.qualifiers:
                    hostTypeDict[isolateHostLookup(seq_feat.qualifiers['isolate'][0])] += 1
                elif 'isolation_source' in seq_feat.qualifiers:
                    hostTypeDict[seq_feat.qualifiers['isolation_source'][0]] += 1
                else: 
                    print("Missing id:", seq_record.id) 
                    notesWithoutHost.append(seq_feat.qualifiers)
                break

    # Verbose print outs and matplotlib:
    # How many was entrex able to query from genbank
    print("Unable to querry the following accession codes: ", set(recordsIds) - set(idsRead), " Successful api querry: {}/{}".format(len(idsRead), len(recordsIds)))
    print(f'\n Host map before further parsing: {hostTypeDict} with size:', sum(x for x in hostTypeDict.values()))
    print(f'\n\n DEBUG: Isolate lookup failures: {isolateLookupFailures} with size:', len(isolateLookupFailures))
    print(f'\n\n DEBUG: Entries with notes without host data: {notesWithoutHost} with size:', len(notesWithoutHost))
    print(f'\n Frequency map before further parsing: {genoTypeDict} with size:', sum(x for x in genoTypeDict.values()))
    plotCounter(genoTypeDict, "pre-processing-genotype-" + (args.p if args.p else "None"), "genotype (raw labels)")
    plotCounter(hostTypeDict, "pre-processing-host-" + (args.p if args.p else "None"), "host (raw labels)")


    ''' --------------- Clean up map for host --------------'''
    species2group = {
            "pig" : "Boar",
            "NoIsolateLookupExists" : "Unknown",
            "swine": "Boar",
            "Homo sapiens": "Human",
            "domestic pig" : "Boar",
            'humna' : 'Human',
            'wild boar' : 'Boar',
            "Rattus norvegicus (wild rat)" : "Rat", 
            "chicken; breeder hen" : "Avian",
            "rabbit" : "Leporidae",
            "Sus scrofa (wild boar)" : "Boar",
            "HepG2/C3A cells" : "Human", # *** These are human cells cancerous test ***
            "Mustela putorius" : "Weasel" ,  # Weasel like?
            "chicken" : "Avian",
            "Rattus tanezumi" : "Rodent",
            "Sus scrofa" : "Boar",
            "Camelus dromedarius" : "Camel",
            "Bos grunniens" : "Bovine",
            "Yunling goat" : "Goat",
            "Macaca mulatta" : "Primate",
            "tree shrew" : "Tree Shrew",
            "cow" : "Bovine",
            "Rhesus monkey serum" : "Primate",
            "Sus scrofa six-month-old" : "Boar",
            "Sus scrofa domesticus" : "Boar",
            "Oryctolagus cuniculus" : "Leporidae" ,
            "wild Boar" : "Boar",
            "Falco tinnunculus" :  "Avian",
            "Rattus norvegicus" : "Rodent",
            "Mustela putorius furo" : "Ferret",
            "Camelus bactrianus" : "Camel",
            "Egretta garzetta (little egret)": "Avian",
            "serum" : "Human", # All serums are human from what I've seen the the database
            "Sus scrofa leucomystax" : "Boar",
            "sparrow" : "Avian",
            "domestic swine" : "Boar",
            "Rattus sp." : "Rodent",
            "Macaca fascicularis" : "Primate",
            "Lepus europaeus (European brown hare)" : "Leporidae",
            "Swine" : "Boar",
            "Tibetan pigs" : "Boar",
            "Microtus arvalis" : "Rodent",
            "Rattus losea" : "Rodent",
            "Homo sapiens PLC/PRF/5 cells" : "Human",
            "Suncus murinus" : "Shrew",
            "swine (domestic pig) serum" : "Boar",
            "human" : "Human",
            "Wild boar" : "Boar",
            "Cervus nippon (wild sika deer)" : "Deer",
            "Herpestes javanicus" : "Mongoose",
            "wolf" : "canine",
            "liver of a pig with progressive weight loss" : "Boar",
            "Oryctolagus cuniculus breed: Rex rabbit": "Leporidae",
            "HepaRG cell culture supernatant" : "Human",
            "Rattus rattus (wild rat)" : "Rodent",
            "Rattus rattus" : "Rodent",
            "human blood plasma" : "Human"
    }

    keys = list(hostTypeDict.keys())
    for i in range(len(hostTypeDict)):
        k = keys[i]
        if k in species2group and species2group[k] != k: # If we know its mapping and its mapping isn't already defined
            hostTypeDict[species2group[k]] += hostTypeDict[k]
            del hostTypeDict[k]
        elif k not in species2group: # unknown entries
            hostTypeDict["unknown"] += hostTypeDict[k]
            del hostTypeDict[k]


    print(f'\n Host map before further parsing: {hostTypeDict} with size:', sum(x for x in hostTypeDict.values()))
    plotCounter(hostTypeDict, "post-processing-host-" + (args.p if args.p else "None"), "host (processed)")

    ''' --------------- Clean up map for genotype --------------'''
    keys = list(genoTypeDict.keys())
    wantedMapping = set(str(x) for x in range(1,8+1))
    wantedMapping.add("None")
    romanNumeralsMap = {"I":"1", "II":"2", "III":"3", "IV":"4", "V":"5", "VI":"6", "VII":"7", "VIII":"8"}
    for i in range(len(genoTypeDict)):
        k = keys[i]

        # If our key is already the one we want, just continue
        if k in wantedMapping:
            continue

        # Move roman numerals I...VIII into 1...8
        flag = 0
        for kR, vR in romanNumeralsMap.items():
            if kR == k: # if the roman number is the key
                genoTypeDict[vR] += genoTypeDict[k] 
                del genoTypeDict[k]
                flag |= 0x1
                break
        if flag:
            continue

        # Move all keys that contain 1,2,3 ... into 1,2,3... respectively
        for j in range(1,8+1): 
            match = str(j)
            if match in k:
                if match != k:
                    genoTypeDict[match] += genoTypeDict[k]
                    del genoTypeDict[k]
                    flag |= 0x01
                break

        # If we couldn't parse it, put it in none (from our ascension example, these misc are just data errors in the database)
        if not flag: 
            genoTypeDict["None"] += genoTypeDict[k]
            del genoTypeDict[k]
            

    print(f'\nFrequency map after mapping {genoTypeDict}', sum(x for x in genoTypeDict.values()))
    print(f'\n\n DEBUG: Entries with notes without genotype data{notesWithoutGenotype} with size:', len(notesWithoutGenotype))
    plotCounter(genoTypeDict, "Processed-genotype-" + (args.p if args.p else "None"), "genotype (processed)")




### -------------------------------------- Main method end --------------------------------------- ##

if __name__ == "__main__":
    args = arg_parse()
    print("Starting program")
    main(args)
