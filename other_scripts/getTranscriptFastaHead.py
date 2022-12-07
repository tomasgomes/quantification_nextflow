from sys import argv

infile = argv[1]
outfile = argv[2]

with open(infile, "r") as inf:
    with open(outfile, "w") as res:
        for line in inf:
            if line.startswith(">"):
                line = line.strip(">").split("|")[-1]
                res.write(">"+line)
            else:
                res.write(line)

