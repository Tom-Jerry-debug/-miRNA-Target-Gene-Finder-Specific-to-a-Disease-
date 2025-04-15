import pandas as pd

def fetch_mirna_targets(mirnas, file_path='data/mirna_target_data.csv'):
    df = pd.read_csv(file_path)
    result = df[df['miRNA'].isin(mirnas)]
    return result