from sys import argv

infile = argv[1]
outfile = argv[2]

with open(infile, "r") as inf:
    with open(outfile, "w") as res:
        for line in inf:
            line = line.split("\t")
            line[-1] = 'transcript_id '+ line[-1].split(";")[0].split(" ")[-1]
            res.write("\t".join(line)+"\n")
