import gseapy as gs
from tableanalyser import get_symbol

def get_ontology_df(topic, background):
    gene_ontology = gs.enrichr(topic, gene_sets='GO_Molecular_Function_2018,GO_Biological_Process_2018,GO_Cellular_Component_2018,Human_Phenotype_Ontology,GTEx_Tissue_Sample_Gene_Expression_Profiles_up,GTEx_Tissue_Sample_Gene_Expression_Profiles_down,Tissue_Protein_Expression_from_Human_Proteome_Map,KEGG_2019_Human,NCI-60_Cancer_Cell_Lines,RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO', cutoff=0.05).results
    return gene_ontology[gene_ontology['Adjusted P-value']<5e-1].loc[:,['Term','Adjusted P-value', 'Gene_set']]

def ensg_to_symbol(serie):
    symbols = []
    for g in serie:
        symbols.append(get_symbol(g))
    return symbols
