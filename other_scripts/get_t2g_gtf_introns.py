from sys import argv
import pandas as pd

gtfin = argv[1]
t2gout = argv[2]
t2g_gn = argv[3]
t2gout_note = argv[4]
intronsout = argv[5]

# parse gtf file
tabnames = pd.read_csv(gtfin, keep_default_na = False, sep = "\t", comment="#")
tabnames = tabnames.loc[tabnames.iloc[:,2]=="exon",:]

gene_transcript = [x[1:] for x in tabnames.iloc[:,8].str.split("; ")]
# get the gene ID
genes = [x[0].split(" ")[1].strip('"') for x in gene_transcript]
# get the transcript ID
transcripts = [x[1].split(' "')[1].strip('"') for x in gene_transcript]
transcripts = [x.split("|")[-1].strip('";') for x in transcripts]
# get the best annotation for the transcript (default is gene name)
annot_gene = [x[1].split(' "')[1].strip('"')for x in gene_transcript]
for i in range(0,len(annot_gene)):
	if "[hs]" in annot_gene[i]:
		annot_gene[i] = annot_gene[i].split(" [hs]")[0].split("|")[-1]
	elif "[nr]" in annot_gene[i]:
		annot_gene[i] = annot_gene[i].split(" [nr]")[0].split("|")[-1]
	elif "|" in annot_gene[i] and not annot_gene[i].startswith("|"):
		annot_gene[i] = annot_gene[i].split("|")[0]
	elif annot_gene[i].startswith("|"):
		annot_gene[i] = annot_gene[i].split("|")[-1].split(".")[0]
	else:
		annot_gene[i] = annot_gene[i].split("|")[-1].split(".")[0]

# get the full annotation for the transcript
annot_full = [x[1].split(' "')[1].strip('"') for x in gene_transcript]

# save it all to a table
t2gnote = pd.DataFrame([transcripts, genes, annot_gene, annot_full]).T

t2gnote.to_csv(t2gout_note, sep = "\t", header=False, index=False)
t2gnote.iloc[:,0:2].to_csv(t2gout, sep = "\t", header=False, index=False)
t2gnote.iloc[:,[0,2]].to_csv(t2g_gn, sep = "\t", header=False, index=False)
t2gnote.iloc[:,0:1].to_csv(intronsout, sep = "\t", header=False, index=False)
