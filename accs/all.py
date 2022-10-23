lines1 = open("ViprBRC-07-04-2022-accession-genes.acc").read().splitlines()
lines2 = open("ViprBRC-07-20-2022-accession-genes.acc").read().splitlines()
lines3 = open("NCBI-07-20-2022-accesssion-genes.acc").read().splitlines()


s = set(lines1 + lines2 + lines3)

with open('all.acc', 'w') as f:
    for line in s:
        f.write("%s\n" % line)


