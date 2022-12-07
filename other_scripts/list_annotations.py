from sys import argv
import pandas as pd

t2gout_note = argv[1]
gene_annot = argv[2]
new_t2g = argv[3]


note_genes = pd.read_csv(t2gout_note, sep = "\t", header=None)
note_genes["annots"] = note_genes.iloc[:,3].str.replace(" [nr]|", "|", regex = False)
note_genes["annots"] = note_genes["annots"].str.replace(" [hs]|", "|", regex = False)
note_genes["annots"] = [x[0] for x in note_genes["annots"].str.split("AMEX60")]

gene_gene = note_genes.iloc[:,[1,2]]
gene_gene["id_name"] = gene_gene.iloc[:,0]+"_"+gene_gene.iloc[:,1]
gene_gene = gene_gene.iloc[:,[0,2]].drop_duplicates(subset=[1])

note_genes = note_genes.iloc[:,[1,4]]
note_genes = note_genes.groupby(1).sum("annots")
note_genes["annots"] = note_genes["annots"].str.strip("|")
note_genes["annots2"] = [",".join(pd.unique(x)) for x in note_genes["annots"].str.split("|")]
note_genes["gene"] = note_genes.index.values

merged_genes = gene_gene.merge(right = note_genes, left_on = 1, right_on = "gene").iloc[:,[0,1,2]]
note_genes = pd.read_csv(t2gout_note, sep = "\t", header=None).iloc[:,0:2]
merged_trans = note_genes.merge(right = merged_genes, left_on = 1, right_on = 1)

merged_trans.to_csv(gene_annot, sep = "\t", header=False, index=False)
merged_trans.iloc[:,[0,2]].to_csv(new_t2g, sep = "\t", header=False, index=False)
