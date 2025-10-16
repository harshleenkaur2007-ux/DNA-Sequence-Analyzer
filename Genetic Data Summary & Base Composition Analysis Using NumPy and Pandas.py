import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt


def load_data(file_path):
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        df = pd.DataFrame(columns=['Sample_ID', 'DNA_Sequence'])
    return df

def add_sequence(df):
    sample_id = input("Enter Sample ID: ")
    dna_seq = input("Enter DNA Sequence: ").upper()
    valid_bases = "ACGT"
    for base in dna_seq:
        if base not in valid_bases:
            print('Invalid DNA Sequence!.Please use only A,C,G,T')
    df = pd.concat([df, pd.DataFrame({'Sample_ID':[sample_id],'DNA_Sequence':[dna_seq]})], ignore_index=True)
    df.to_csv(file_path, index=False)
    print(f"Added sequence for {sample_id} and saved to {file_path}")
    return df

def analyze_sequences(df):
    A,C,T,G,GC = [],[],[],[],[]
    for seq in df['DNA_Sequence']:
        A.append(seq.count('A'))
        C.append(seq.count('C'))
        T.append(seq.count('T'))
        G.append(seq.count('G'))
        GC.append(round(((seq.count('G') + seq.count('C')) / len(seq)) * 100, 2))
    df['A_Count'] = A; df['C_Count'] = C; df['T_Count'] = T; df['G_Count'] = G; df['GC_Content(%)'] = GC
    df['Category'] = ['High GC' if g>60 else 'Low GC' if g<40 else 'Normal' for g in GC]
    print('Analysis Done')
    return df

def cluster_sequences(df, n_clusters=3):
    if len(df) < n_clusters:
        df['Cluster'] = -1
        print("Not enough sequences for clustering.")
        return df

    X = df[['A_Count','T_Count','G_Count','C_Count','GC_Content(%)']]
    df['Cluster'] = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(X)

    print(f"Clustered into {n_clusters} groups.\n")
    plt.scatter(df['G_Count'], df['C_Count'], c=df['Cluster'], cmap='viridis', s=100)
    plt.xlabel('G Count'); plt.ylabel('C Count'); plt.title('DNA Sequence Clusters')
    plt.colorbar(label='Cluster'); plt.show()
    print(df[['Sample_ID','DNA_Sequence','Cluster']])
    return df

def display_results(df):
    if df.empty:
        print("No data.")
        return
    print(df)
    print(f"Average GC Content: {np.mean(df['GC_Content(%)']):.2f}%")


print("Welcome to the Simple DNA Sequence Analyzer!")
print("This tool lets you add DNA sequences, analyze base counts and GC content,")
print("cluster sequences into groups, and visualize or save the results.\n")

file_path = 'dna_data.csv'
output_file = 'dna_results.csv'
df = load_data(file_path)
    
while True:
    print("\n===== DNA Sequence Analyzer Menu =====")
    print("1. Add a new DNA sequence")
    print("2. Analyze sequences (base counts + GC content)")
    print("3. Cluster sequences")
    print("4. Display results")
    print("5. Exit")
    choice = input("Enter your choice (1-6): ")
    if choice=='1':
        df = add_sequence(df)
    elif choice=='2':
        df = analyze_sequences(df)
    elif choice=='3':
        if 'GC_Content(%)' not in df or df['GC_Content(%)'].isnull().any():
            print("Running analysis first...")
            df = analyze_sequences(df)
        n = input("Enter number of clusters (default 3): ")
        n = int(n) if n.isdigit() else 3
        df = cluster_sequences(df, n)
    elif choice=='4':
        display_results(df)
    elif choice=='5':
        print("Exiting program. Goodbye!")
        break
    else:
        print("Invalid choice! Please enter a number between 1 and 6.")













