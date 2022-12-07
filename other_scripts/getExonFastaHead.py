from sys import argv
import pandas as pd

infile = argv[1]
outfile = argv[2]
exong = argv[3]

tr_dic = {}
genel = []
exonl = []
with open(infile, "r") as inf:
    with open(outfile, "w") as res:
        for line in inf:
            if line.startswith(">"):
                line = line.strip(">").split("|")[-1].strip("\n")
                if not line in tr_dic.keys():
                        tr_dic[line] = 1
                else:
                        tr_dic[line]+=1
                line+="."+str(tr_dic[line])
                genel.append(line.split(".")[0])
                exonl.append(line)
                res.write(">"+line+"\n")
            else:
                res.write(line)

pd.DataFrame(list(zip(exonl, genel))).to_csv(exong, sep="\t", header=False, index=False)
