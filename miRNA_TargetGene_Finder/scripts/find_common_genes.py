import pandas as pd

def find_common_genes(mirna_targets, disease, file_path='data/disease_gene_data.csv'):
    disease_df = pd.read_csv(file_path)
    disease_genes = disease_df[disease_df['Disease'] == disease]['Gene_Symbol'].tolist()
    
    common = mirna_targets[mirna_targets['Target_Gene'].isin(disease_genes)]
    return common