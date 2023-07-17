#A Colab file written in Python to create a graph plotting the UMI expressions of genes on a given chromosome (in this case the Y) in order of chromosonal location

#importing necessary packages
import pandas as pd
import matplotlib.pyplot as plt

#reading the files, replace with those specific for your project
genes = pd.read_csv("geneposition.csv")
umis = pd.read_csv("Names and Total UMI for Y Chromsome Genes.csv")

#removing genes not associated with hgnc symbols for Y chromosome-linked genes
rows_to_drop = []
for row in genes.itertuples():
    if pd.isna(row.hgnc_symbol):
        rows_to_drop.append(row.Index)
genes = genes.drop(rows_to_drop, axis=0)

#renaming column for consistency
genes = genes.rename(columns={'hgnc_symbol': 'Gene'})

#dropping duplicate entries
genes.drop_duplicates(subset=['Gene'], inplace = True)
umis.drop_duplicates(subset=['Gene'], inplace = True)

#merging the two files to create an integrated one
merged_df = pd.merge(genes, umis, on='Gene',how='inner')

#dropping unneeded columns
m = merged_df.drop(['ensembl_gene_id','UMI','Unnamed: 0'],axis=1)

#sort in order of chromosonal position
m.sort_values('start_position', ascending=True,inplace=True)

#remove genes that have 0 UMI (hard to display on graph)
m = m[m['Total UMI'] != 0]

#scaling data to fit better on graph
m_root = m.copy()
m_root['start_position'] = m['start_position'] * 0.000001

#creating graph (change these in accordance with specific needs of project)
pos = m_root['start_position']
umi = m_root['Total UMI']
plt.bar(pos, umi, width = 0.1)
plt.xlabel("Relative starting position of gene")
plt.ylabel("Y UMI (normalized)")
plt.title("Y UMI across genes located at different starting positions along the Y-chromosome")
plt.show()



