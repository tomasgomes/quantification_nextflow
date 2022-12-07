# faidx --split-files ../AmexG_v6.DD.corrected.round2.chr.fa
# faidx ../AmexG_v6.DD.corrected.round2.chr.fa -i chromsizes > AmexG_v6.DD.corrected.round2.chr.sizes

# python cutChroms.py AmexG_v6_chr/AmexG_v6.DD.corrected.round2.chr.sizes ../AmexT_v47.FULL_corr_chr.gtf AmexG_v6_chr/AmexT_v47.FULL_corr_chr_cut.gtf AmexG_v6_chr/cutcoord.txt AmexG_v6_chr
# cat *.fa > AmexG_v6.chr_cut.fa


from sys import argv
from os import system

chrsizes = argv[1] # chromosome sizes
gtffile = argv[2] # input gtf for gene coordinates
outgtf = argv[3] # new gene coordinates
outcut = argv[4] # cut coordinates
outdir = argv[5] # output directory for the new cut chromosomes
maxSize = int(argv[6]) # 500000000
spSize = int(argv[7]) # 2000
cutdic = {}

# read chr file into dic
sizesdic = {}
with open(chrsizes, "r") as f:
  for line in f:
    line = line.split("\t")
    sizesdic[line[0]] = int(line[1])


# read genes end coord in GTF into dictionary per chromosome ({"chr1": {"geneA": 100}})
genedic = {}
# also read all GTF entries into a dictionary, splitting the lines ({"chr1": {"geneA": [...], "exon1A": [...]}})
gtfdic = {}
with open(gtffile, "r") as f:
  for line in f:
    line = line.split("\t")
    line[3] = int(line[3])
    line[4] = int(line[4])
    linegeneid = line[8].split(";")[0].split(" ")[1].strip('"')

    if line[2]=="gene":
      if line[0] in genedic.keys():
        genedic[line[0]][linegeneid] = [line[3], line[4]]
      else:
        genedic[line[0]] = {linegeneid: [line[3],line[4]]}
        
    if line[0] in gtfdic.keys():
      if linegeneid in gtfdic[line[0]].keys():
        gtfdic[line[0]][linegeneid].append(line)
      else:
        gtfdic[line[0]][linegeneid] = [line]
    else:
      gtfdic[line[0]] = {linegeneid: [line]}
        

# for each chr (based on size dic)
for chrom in sizesdic.keys():
  print(chrom)
  it = 1
  en = 0
  s = sizesdic[chrom]
  ## check size. if big: (while chrsize>maxSize:)
  ### find gene in that chromosome with end closest to but within first break
  ### subtract maxSize to the chr size, all gene coordinates (gene dic) or the original chr, and all entries for that chr in the GTF dic
  ### put all genes up to that in a new chromosome (e.g. chr1p_1) in the gene dic, remove from old
  ### put all GTF records up to that in a new chromosome (e.g. chr1p_1) in the GTF dic, remove from old
  ### add coordinates to cutdic (used to cut the fasta file)
  while s>maxSize:
    newchr = chrom+"s"+str(it)
    print(newchr)
    
    wthg = {g:genedic[chrom][g] for g in genedic[chrom].keys() if (genedic[chrom][g][1]-(en+maxSize))<0}
    genedic[newchr] = wthg
    
    cutg = [genedic[chrom][g][0] for g in genedic[chrom].keys() if (genedic[chrom][g][0]-(en+maxSize))<0 and (genedic[chrom][g][1]-(en+maxSize))>0]
    xxx = [g for g in genedic[chrom].keys() if (genedic[chrom][g][0]-(en+maxSize))<0 and (genedic[chrom][g][1]-(en+maxSize))>0]
    
    st = en+1
    if len(cutg)>0:
      while len(cutg)>0:
        en = min(cutg)-spSize
        cutg = [genedic[chrom][g][0] for g in genedic[chrom].keys() if (genedic[chrom][g][0]-(en))<0 and (genedic[chrom][g][1]-(en))>0]
        badg = [g for g in genedic[chrom].keys() if (genedic[chrom][g][0]-(en))<0 and (genedic[chrom][g][1]-(en))>0]
        if len(badg)>0:
          for gg in badg:
            genedic[newchr].pop(gg)
    else:
      en = st-1+maxSize
    cutdic[newchr] = [str(st), str(en)]
    
    for g in wthg.keys(): # remove those genes
      genedic[chrom].pop(g)
    
    gtfdic[newchr] = {g:gtfdic[chrom][g] for g in wthg.keys()}
    for g in gtfdic[newchr].keys():
      for i in range(0, len(gtfdic[newchr][g])):
        gtfdic[newchr][g][i][0] = newchr
        gtfdic[newchr][g][i][3] = str(gtfdic[newchr][g][i][3]-st+1)
        gtfdic[newchr][g][i][4] = str(gtfdic[newchr][g][i][4]-st+1)
    
    # update condition and it
    s = s-en+st
    print(s)
    it+=1
      
  ## else:
  ### if it!=1
  #### add coordinates to cutdic (used to cut the fasta file)
  #### put all GTF records up to that in a new chromosome (e.g. chr1p_1) in the GTF dic, remove from old
  else:
    if it!=1:
      newchr = chrom+"s"+str(it)
      genedic[newchr] = genedic[chrom]
      
      st = en+1
      cutdic[newchr] = [str(st), str(sizesdic[chrom])]
      
      gtfs = {g:gtfdic[chrom][g] for g in genedic[newchr].keys()}
      gtfdic[newchr] = gtfs
      for g in gtfdic[newchr].keys():
        for i in range(0, len(gtfdic[newchr][g])):
          gtfdic[newchr][g][i][0] = newchr
          gtfdic[newchr][g][i][3] = str(gtfdic[newchr][g][i][3]-st+1)
          gtfdic[newchr][g][i][4] = str(gtfdic[newchr][g][i][4]-st+1)
      gtfdic.pop(chrom)
      
    else:
      cutdic[chrom] = [str(1), str(s)]
      for g in gtfdic[chrom].keys():
        for i in range(0, len(gtfdic[chrom][g])):
          gtfdic[chrom][g][i][3] = str(gtfdic[chrom][g][i][3])
          gtfdic[chrom][g][i][4] = str(gtfdic[chrom][g][i][4])


print("Creating new GTF file...")
# write gtf
with open(outgtf,"w") as resgtf:
  for c in gtfdic.keys():
    for g in gtfdic[c].keys():
      for l in gtfdic[c][g]:
        resgtf.write("\t".join(l))


print("Creating file for fasta cuts...")
# write cut coord
with open(outcut,"w") as rescut:
  for c in cutdic.keys():
    rescut.write(c.split("s")[0]+":"+"-".join(cutdic[c])+"\n")


print("Cutting fasta...")
for c in cutdic.keys():
  print(c, end = " ")
  coord = c.split("s")[0]+":"+"-".join(cutdic[c])
  system("samtools faidx AmexG_v6.DD.corrected.round2.chr.fa "+coord+" > "+outdir+"/AmexG_v6."+c+".fa")
  system("sed -i 's/"+coord+"/"+c+"/g' "+outdir+"/AmexG_v6."+c+".fa")

system("cat "+outdir+"/*.fa > "+outdir+"/AmexG_v6.chr_cut.fa")




