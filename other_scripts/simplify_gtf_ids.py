from sys import argv

infile = argv[1]
outfile = argv[2]

with open(outfile, "w") as res:
    with open(infile, "r") as inp:
        for line in inp:
            line = line.split("\t")
            line[8] = line[8].split("; ")
            line[8] = "; ".join(line[8][0:3])
            res.write("\t".join(line)+"\n")
