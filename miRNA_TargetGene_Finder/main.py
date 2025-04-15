from scripts.fetch_genes import fetch_mirna_targets
from scripts.find_common_genes import find_common_genes
from scripts.fetch_sequences import fetch_gene_description
import pandas as pd

def main():
    mirnas = ['miR-1', 'miR-133']
    disease = 'cardiac hypertrophy'

    targets = fetch_mirna_targets(mirnas)
    print("\nTarget Genes for input miRNAs:\n", targets)

    common_genes = find_common_genes(targets, disease)
    print("\nCommon Genes with", disease, ":\n", common_genes)

    # Fetch descriptions
    descriptions = []
    for gene in common_genes['Target_Gene']:
        desc = fetch_gene_description(gene)
        descriptions.append(desc)

    common_genes['Gene_Description'] = descriptions
    print("\nFinal Annotated Results:\n", common_genes)

    # Save result
    common_genes.to_csv('output/result.csv', index=False)

if __name__ == "__main__":
    main()