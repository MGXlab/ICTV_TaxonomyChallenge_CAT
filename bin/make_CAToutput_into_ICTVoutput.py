# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 10:30:05 2024

@author: Ernestina Hauptfeld
"""

ICTV_ranks=['clade', 'subrealm', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 
            'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 
            'genus', 'subgenus', 'species']


def make_ncbi2ictv_dict(input_file):
    ncbi2ictv={}
    
    tmp=open(input_file).read().strip().split('\n')[3:]
        
    for acc in tmp:
        acc=acc.split('\t')
        
        taxid=acc[2]
        accession=acc[0]
        ictv_species=acc[4].split(';')
        ncbi_species=acc[6].split(';')
        ncbi_ranks=acc[7].split(';')
        ncbi_taxids=acc[5].split(';')
        
        if taxid not in ncbi2ictv:
            ncbi2ictv[taxid]={'accessions': [],
                              'ictv': [],
                              'ncbi': [],
                              'ncbi_ranks': [],
                              'ncbi_taxids': []}
        
        ncbi2ictv[taxid]['accessions'].append(accession)
        ncbi2ictv[taxid]['ictv'].append(ictv_species)
        ncbi2ictv[taxid]['ncbi'].append(ncbi_species)
        ncbi2ictv[taxid]['ncbi_ranks'].append(ncbi_ranks)
        ncbi2ictv[taxid]['ncbi_taxids'].append(ncbi_taxids)
        
        for rank in ncbi_ranks:
            if rank in ICTV_ranks:
                short_lineage=ncbi_taxids[:ncbi_ranks.index(rank)+1]
                new_taxid=short_lineage[-1]
                short_names=ncbi_species[:ncbi_ranks.index(rank)+1]
                short_ranks=ncbi_ranks[:ncbi_ranks.index(rank)+1]
                short_ictv=ictv_species[:ICTV_ranks.index(rank)+1]
                if new_taxid not in ncbi2ictv:
                    ncbi2ictv[new_taxid]={'accessions': [],
                                      'ictv': [],
                                      'ncbi': [],
                                      'ncbi_ranks': [],
                                      'ncbi_taxids': []}
                ncbi2ictv[new_taxid]['ictv'].append(short_ictv)
                ncbi2ictv[new_taxid]['ncbi'].append(short_names)
                ncbi2ictv[new_taxid]['ncbi_ranks'].append(short_ranks)
                ncbi2ictv[new_taxid]['ncbi_taxids'].append(short_lineage)
        
    return ncbi2ictv


def make_c2c_dict(c2c_file):
    c2c={}
    
    tmp=open(c2c_file).read().strip().split('\n')[1:]
    
    for contig in tmp:
        contig=contig.replace('*','').split('\t')
        contig_id=contig[0]
        c2c[contig_id]={'lineage': [],
                     'ranks': [],
                     'scores': [],
                     'names': []}
        
        if contig[1]!='no taxid assigned':
            lineage=contig[3].split(';')
            scores=contig[4].split(';')
            
            names_ranks=contig[5:]
            names=[]
            ranks=[]
            for r in names_ranks:
                name=r.split(' (')[0]
                rank=r.split(' (')[1].split('):')[0]
                names.append(name)
                ranks.append(rank)
            
            c2c[contig_id]['lineage']=lineage
            c2c[contig_id]['scores']=scores
            c2c[contig_id]['ranks']=ranks
            c2c[contig_id]['names']=names
            
    return c2c
        
def find_LCA(list_of_lineages):
    lca=[]
    overlap = set.intersection(*map(set, list_of_lineages))

    for taxid in list_of_lineages[0]:
        if taxid in overlap:
            lca.append(taxid)
        else:
            break
    return lca


def make_challenge_output_direct(c2c_dict, output_file):
    output_dict={}
    for contig in c2c_dict:
        output_dict[contig]={'names': len(ICTV_ranks)*['NA'],
                             'scores': len(ICTV_ranks)*['NA']}
        for rank in ICTV_ranks:
            if rank in c2c_dict[contig]['ranks']:
                c2c_index=c2c_dict[contig]['ranks'].index(rank)
                output_index=ICTV_ranks.index(rank)
                
                output_dict[contig]['names'][output_index]=c2c_dict[contig]['names'][c2c_index]
                output_dict[contig]['scores'][output_index]=c2c_dict[contig]['scores'][c2c_index]
                
    
    with open(output_file, 'w') as outf:
        outf.write('SequenceID,realm,realm_score,subrealm,subrealm_score')
        for rank in ICTV_ranks[2:]:
            outf.write(f',{rank},{rank}_score')
        outf.write('\n')
        
        for contig in output_dict:
            outf.write(contig)
            for i in range(0,15):
                outf.write(f",{output_dict[contig]['names'][i]},{output_dict[contig]['scores'][i]}")
            outf.write('\n')
        
    return


def make_challenge_output_mapped(c2c_dict, ncbi2ictv,output_file):
    n=0
    output_dict={}
    for contig in c2c_dict:
        n+=1
        if n%1000==0:
            print(f'done with {n} contigs!')
        output_dict[contig]={'names': len(ICTV_ranks)*['NA'],
                             'scores': len(ICTV_ranks)*['NA']}
        mapped=False
        for i in range(1,len(c2c_dict[contig]['ranks'])+1):
            if not mapped and c2c_dict[contig]['ranks'][-i] in ICTV_ranks:
                mapped=True
                taxid=c2c_dict[contig]['lineage'][-i]
                
                if taxid in ncbi2ictv:
                    ictv_lineage=find_LCA(ncbi2ictv[taxid]['ictv'])
                    ncbi_lineage=find_LCA(ncbi2ictv[taxid]['ncbi'])
                    ncbi_tmp=find_LCA(ncbi2ictv[taxid]['ncbi_ranks'])
                    ncbi_ranks=ncbi_tmp[:len(ncbi_lineage)+1]
                    
                    for rank in ictv_lineage:
                        ICTV_rank=ICTV_ranks[ictv_lineage.index(rank)]
                        if ICTV_rank in ncbi_ranks:
                            output_dict[contig]['names'][ICTV_ranks.index(ICTV_rank)]=ictv_lineage[ICTV_ranks.index(ICTV_rank)]
                           
    
    print('Writing output file...')
    with open(output_file, 'w') as outf:
        outf.write('SequenceID,realm,realm_score,subrealm,subrealm_score')
        for rank in ICTV_ranks[2:]:
            outf.write(f',{rank},{rank}_score')
        outf.write('\n')
        
        for contig in output_dict:
            outf.write(contig)
            for i in range(0,15):
                outf.write(f",{output_dict[contig]['names'][i]},{output_dict[contig]['scores'][i]}")
            outf.write('\n')
        
        
        
    return


if __name__=='__main__':    
    ncbi2ictv=make_ncbi2ictv_dict('../ICTV39_NCBI202412_per_accession.tsv')
    c2c=make_c2c_dict('../20241217_cat.c2c.names.txt')
    make_challenge_output_direct(c2c, '../results/CAT_pack_ICTV_challenge_submission1.csv')
    make_challenge_output_mapped(c2c, ncbi2ictv, '../results/CAT_pack_ICTV_challenge_submission2.csv')