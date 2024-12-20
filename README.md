# CAT_pack participation for ICTV Taxonomy Challenge

This is the git repository containing the results of CAT_pack v6.0.1 annotating the contigs of the ICTV taxonomy challenge.

## Contributors
Ernestina Hauptfeld: CAT_pack development/maintenance, performed this analysis  
F.A. Bastiaan von Meijenfeldt: CAT_pack development/maintenance  
Bas E. Dutilh: Supervision  

## CAT_pack
The contig annotation tool CAT (i) predicts proteins on contigs using prodigal, (ii) runs a diamond search of the predicted proteins against the NR database, (iii) assigns a taxonomic lineage to each predicted protein using the last common ancestor of all diamond hits within 10% of the top bitscore, and (iv) assigns a taxonomic lineage to each contig based on the lineages of all ORFs on the contig.
[Link to the CAT paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1817-x)
[Link to the CAT github](https://github.com/MGXlab/CAT_pack)

We ran CAT_pack version 6.0.1, using diamond

CAT_pack contigs -c /home/tina/Documents/ICTV_challenge/dataset_challenge/ICTV_all_contigs.fasta -t CAT_prepare/20241212_CAT_nr/tax/ -d CAT_prepare/20241212_CAT_nr/db/ -n 96 -o ICTV/cat_with_ncbi_uncorrected/20241217_cat_ncbi_uncorrected --path_to_diamond CAT_prepare/20241212_CAT_nr/diamond
