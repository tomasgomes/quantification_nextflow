import sys

chrsizes = sys.argv[1] # chrName lenght
coord = sys.argv[2]

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
chrname = coord.split(":")[0]
coord1 = coord.split(":")[1].split("-")[0]
coord2 = coord.split(":")[1].split("-")[1]
if "s" in chrname:
  coord1 = str((size_factor_dic[chrname] + int(coord1)))
  coord2 = str((size_factor_dic[chrname] + int(coord2)))
  chrname = chrname[:len(chrname)-2]
else:
  coord1 = str((int(coord1)))
  coord2 = str((int(coord2)))

print(chrname+":"+str(coord1)+"-"+str(coord2))
