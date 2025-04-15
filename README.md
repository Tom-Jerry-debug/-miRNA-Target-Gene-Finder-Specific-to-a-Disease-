# -miRNA-Target-Gene-Finder-Specific-to-a-Disease-
✅ Step 1: Set Up Your Environment
🔧 Install Required Libraries
Open your terminal or Anaconda Prompt and install:

bash
Copy
Edit
pip install biopython pandas requests matplotlib
If you’re using Jupyter Notebook:

python
Copy
Edit
!pip install biopython pandas requests matplotlib

✅ Step 2: Set Up Directory Structure
Create a folder like this:

css
Copy
Edit
miRNA_TargetGene_Finder/
│
├── data/
│   ├── mirna_target_data.csv
│   └── disease_gene_data.csv
│
├── output/
│   └── result.csv
│
├── scripts/
│   ├── fetch_genes.py
│   ├── find_common_genes.py
│   ├── fetch_sequences.py
│
├── main.py
└── requirements.txt
You can manually create this folder, or use Python:

python
Copy
Edit
import os

folders = ['data', 'output', 'scripts']
for folder in folders:
    os.makedirs(folder, exist_ok=True)
    
✅ Step 3: Collect Sample Data
📑 mirna_target_data.csv

miRNA	Target_Gene
miR-1	IGF1
miR-1	MYH7
miR-133	MYH7
miR-133	ACTA1
📑 disease_gene_data.csv

Disease	Gene_Symbol
cardiac hypertrophy	IGF1
cardiac hypertrophy	MYH7
cardiac hypertrophy	TNNT2
👉 We’ll use this sample data to test the script before going online.


✅ Step 4: Set Up Entrez (NCBI API) in Biopython
You need an NCBI API key for large queries. Create one at https://www.ncbi.nlm.nih.gov/account/

Then configure Entrez:

python
Copy
Edit
from Bio import Entrez

Entrez.email = "your_email@example.com"
Entrez.api_key = "YOUR_API_KEY"

✅ Step 5: Create the Core Python Scripts
📜 scripts/fetch_genes.py
Fetch target genes for given miRNAs:

python
Copy
Edit
import pandas as pd

def fetch_mirna_targets(mirnas, file_path='data/mirna_target_data.csv'):
    df = pd.read_csv(file_path)
    result = df[df['miRNA'].isin(mirnas)]
    return result
📜 scripts/find_common_genes.py
Find common genes between miRNA targets and disease-associated genes:

python
Copy
Edit
import pandas as pd

def find_common_genes(mirna_targets, disease, file_path='data/disease_gene_data.csv'):
    disease_df = pd.read_csv(file_path)
    disease_genes = disease_df[disease_df['Disease'] == disease]['Gene_Symbol'].tolist()
    
    common = mirna_targets[mirna_targets['Target_Gene'].isin(disease_genes)]
    return common
📜 scripts/fetch_sequences.py
Fetch gene info from NCBI:

python
Copy
Edit
from Bio import Entrez

Entrez.email = "your_email@example.com"
Entrez.api_key = "YOUR_API_KEY"

def fetch_gene_description(gene_symbol):
    handle = Entrez.esearch(db="gene", term=gene_symbol + "[sym]")
    record = Entrez.read(handle)
    handle.close()
    
    if record['IdList']:
        gene_id = record['IdList'][0]
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        records = Entrez.read(handle)
        handle.close()
        description = records[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
        return description
    else:
        return "Description not found"
        
✅ Step 6: Main Controller Script
📜 main.py
This glues everything together.

python
Copy
Edit
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
    
✅ Step 7: Run the Program
In terminal:

bash
Copy
Edit
python main.py

✅ It will:

Fetch target genes for miRNAs

Find common genes associated with a disease

Fetch descriptions from NCBI

Save the final results in output/result.csv

✅ Step 8: (Optional) Visualize Results
Using a simple bar graph:

python
Copy
Edit
import matplotlib.pyplot as plt

df = pd.read_csv('output/result.csv')
df.groupby('miRNA')['Target_Gene'].count().plot(kind='bar')
plt.title('Number of Common Genes per miRNA')
plt.ylabel('Count')
plt.show()

📌 Summary:
✅ Created a clean, modular, project
✅ Integrated Biopython’s Entrez for live gene info
✅ Built easy-to-edit CSV-based datasets
✅ Results saved in structured CSV
✅ Visualized data
✅ 100% expandable for larger datasets and diseases



