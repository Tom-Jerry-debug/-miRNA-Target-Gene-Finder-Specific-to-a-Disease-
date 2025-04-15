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