#/usr/bin/env sh

echo VIPRBRC is Dataset 1, NCBI is dataset 2
# ./compareAccs.py -p1 ./accs/ViprBRC-07-04-2022-accession-genes.acc -p2 ./accs/NCBI-07-20-2022-accesssion-genes.acc | tee logs/compare_ViprBRC_NCBI.log
./compareAccs.py -p1 ./accs/ViprBRC-07-20-2022-accession-genes.acc -p2 ./accs/NCBI-07-20-2022-accesssion-genes.acc | tee logs/compare_ViprBRC_NCBI.log