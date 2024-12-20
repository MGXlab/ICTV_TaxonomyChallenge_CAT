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

We ran CAT_pack version 6.0.1, which uses prodigal version 2.6.3 and diamond version 2.1.10.164.  
We compared to an NR database, downloaded on 12th of December, 2024. The ICTV taxonomy used is consistent with ICTV MSL39.  

```
# Running CAT
CAT_pack contigs -c ICTV_all_contigs.fasta -t ./20241212_CAT_nr_website/tax/ -d ./20241212_CAT_nr_website/db/ -n 96 -o ICTV/cat_output/20241217_raw_CAT_results --path_to_diamond ./20241212_CAT_nr_website/diamond
python3 bin/make_CAToutput_into_ICTVoutput.py
```

## Results

*CAT_pack_ICTV_challenge_submission1.csv*: Raw CAT_pack output reformatted into ICTV challenge format.    

*CAT_pack_ICTV_challenge_submission2.csv*: CAT_pack output mapped to ICTV taxonomy via the mapping file *ICTV39_NCBI202412_per_accession.tsv*.  

The mapping file *ICTV39_NCBI202412_per_accession.tsv* was created combining information from the [ICTV Virus Metadata Resource MSL39 v4](https://ictv.global/sites/default/files/VMR/VMR_MSL39.v4_20241106.xlsx) and the taxonomy files of the NR database on the 12th of December. 
