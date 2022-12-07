import sys

chrsizes = sys.argv[1] # chrName lenght
bed = sys.argv[2]
output = sys.argv[3]

# open new chr sizes and compile size to add to each segment
size_factor_dic = {}
with open(chrsizes, "r") as f:
  for line in f:
    line = line.strip("\n").split("\t")
    if "s" in line[0]:
      if line[0].endswith("s1"):
        size_factor_dic[line[0]] = 0
        size_factor_dic[line[0].replace("s1", "s2")] = int(line[1])
      else:
        last_s = line[0][len(line[0])-1:(len(line[0]))]
        last = str(int(last_s)+1)
        size_factor_dic[line[0].replace("s"+last_s, "s"+last)] = int(line[1]) + size_factor_dic[line[0]]
    else:
      size_factor_dic[line[0]] = 0

# change chromosome names and sizes
with open(output, "w") as res:
  with open(bed, "r") as f:
    for line in f:
      line = line.split("\t")
      if "s" in line[0]:
        line[1] = str((size_factor_dic[line[0]] + int(line[1]))+1)
        line[2] = str((size_factor_dic[line[0]] + int(line[2]))+1)
        line[0] = line[0][:len(line[0])-2]
      else:
        line[1] = str((int(line[1])+1))
        line[2] = str((int(line[2])+1))
      
      newline = "\t".join(line)
      res.write(newline)
