genotypeToHosts = {
    "HEV1": {"Human"},
    "HEV2": {"Human"},
    "HEV3": {"Human", "Boar", "Deer", "Primate", "Mongoose", "Leporidae"},
    "HEV3ra": {"Human", "Boar", "Deer", "Primate", "Mongoose", "Leporidae"},
    "HEV4": {"Human", "Boar", "Primate"},
    "HEV5": {"Boar"},
    "HEV6": {"Boar"},
    "HEV7": {'Camel'},
    "HEV8": {'Camel'},
}

''' --------------- Clean up map for host --------------'''
species2group = {
    "Unknown": "Unknown",
    "Avian": "Avian",
    "Swine": "Boar",
    "Boar": "Boar",
    "pig": "Boar",
    "Bat": "Bat",
    "Leporidae": "Leporidae",
    "NoIsolateLookupExists": "Unknown",
    "swine": "Boar",
    "Homo sapiens": "Human",
    "domestic pig": "Boar",
    'humna': 'Human',
    "Human": "Human",
    'wild boar': 'Boar',
    "Rattus norvegicus (wild rat)": "Rodent",
    "chicken; breeder hen": "Avian",
    "rabbit": "Leporidae",
    "Sus scrofa (wild boar)": "Boar",
    "HepG2/C3A cells": "Human",  # *** These are human cells cancerous test ***
    "Mustela putorius": "Weasel",  # Weasel like?
    "chicken": "Avian",
    "Rattus tanezumi": "Rodent",
    "Sus scrofa": "Boar",
    "Camelus dromedarius": "Camel",
    "Bos grunniens": "Bovine",
    "Yunling goat": "Goat",
    "Macaca mulatta": "Primate",
    "tree shrew": "Shrew",
    "cow": "Bovine",
    "Rhesus monkey serum": "Primate",
    "Sus scrofa six-month-old": "Boar",
    "Sus scrofa domesticus": "Boar",
    "Oryctolagus cuniculus": "Leporidae",
    "wild Boar": "Boar",
    "Falco tinnunculus": "Avian",
    "Rattus norvegicus": "Rodent",
    "Mustela putorius furo": "Ferret",
    "Camelus bactrianus": "Camel",
    "Egretta garzetta (little egret)": "Avian",
    "serum": "Human",  # All serums are human from what I've seen the database
    "Sus scrofa leucomystax": "Boar",
    "sparrow": "Avian",
    "domestic swine": "Boar",
    "Rattus sp.": "Rodent",
    "Macaca fascicularis": "Primate",
    "Lepus europaeus (European brown hare)": "Leporidae",
    "Tibetan pigs": "Boar",
    "Microtus arvalis": "Rodent",
    "Rattus losea": "Rodent",
    "Homo sapiens PLC/PRF/5 cells": "Human",
    "Suncus murinus": "Shrew",
    "swine (domestic pig) serum": "Boar",
    "human": "Human",
    "Wild boar": "Boar",
    "Cervus nippon (wild sika deer)": "Deer",
    "Herpestes javanicus": "Mongoose",
    "wolf": "Canine",
    "liver of a pig with progressive weight loss": "Boar",
    "Oryctolagus cuniculus breed: Rex rabbit": "Leporidae",
    "HepaRG cell culture supernatant": "Human",
    "Rattus rattus (wild rat)": "Rodent",
    "Rattus rattus": "Rodent",
    "human blood plasma": "Human",
    "stools of infected French soldiers from 1983 outbreak": "Human",
    "domestic pig serum": "Boar",
    "fecal sample of experimentally infected swine": "Boar",
    "Serum": "Human",
    "culture supernatant": "Unknown",
    "human serum": "Human",
    "wild deer": "Deer",
    "patient serum": "Human",
    "a patient with hepatitis E": "Human",
    "wild boar liver": "Boar",
    "wild mongoose serum": "Mongoose",
    "swine serum": "Boar",
    "feces": "Human",
    "sausage": "Boar",
    "wild boar hybrid": "Boar",
    "blood": "Unknown",
    "liver": "Unknown",
    "Eptesicus serotinus": "Bat",
    "Rattus norvegicus (Norway rat)": "Rodent",
    "Ochotona curzoniae and Apodemus peninsulae": "Hare",
    "Apodemus chevrieri": "Rodent",
    "Eothenomys melanogaster": "Rodent",
    "Gallus gallus": "Avian",
    "Chicken": "Avian",
    "silkie": "Avian"
}

manualMap = {"KX513953": "Bat",  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2168470/
             "LY647099": "Leporidae",  # https://www.ncbi.nlm.nih.gov/nuccore/LY647099
             "AB425831": "Human",
             # https://www.ncbi.nlm.nih.gov/nuccore/AB425831 https://pubmed.ncbi.nlm.nih.gov/6889516/
             "AB480825": "Human",  # https://pubmed.ncbi.nlm.nih.gov/19369433/
             "LC061267": "Human",  # PCL/PRF/5 cell line https://www.ncbi.nlm.nih.gov/nuccore/LC061267
             "AB362839": "Human",  # Same as above https://www.ncbi.nlm.nih.gov/nuccore/AB362839
             "AB425830": "Human",  # Same as above https://pubmed.ncbi.nlm.nih.gov/18620009/
             "AB362841": "Human",  # Same as above https://pubmed.ncbi.nlm.nih.gov/18620009/
             "AB362840": "Human",  # A549 cells https://www.ncbi.nlm.nih.gov/nuccore/AB362840
             "AJ272108": "Human",
             # likely human from patient keyword # https://www.microbiologyresearch.org/content/journal/jgv/10.1099
             # /0022-1317-81-7-1675#tab2
             }
