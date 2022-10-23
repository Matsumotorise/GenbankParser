#/usr/src/env sh
#./query_hepA_nuc.py -p ./accs/ViprBRC-07-04-2022-accession-genes.acc | tee logs/VIPRBRC.txt
#./query_hepA.py -p ./accs/sequences.acc tee | tee logs/ncbi_output.log

./query_hepA_nuc.py -p ./accs/ViprBRC-07-04-2022-accession-genes.acc | tee logs/VIPBRC-NCBI.txt

#./accToFasta.py -p ./accs/VIPRBRR_NCBI-accession-genes-07-20-2022.acc

