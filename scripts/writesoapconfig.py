import sys

#USAGE python3 writesoapconfig.py <outfilename> <forward_readfile> <reverse_readfile>

outfile_name = sys.argv[1]
fReadname = sys.argv[2]
rReadname = sys.argv[3]

with open("metaAndconfig/config.yaml") as config, open(outfile_name, "w") as outfile:
    read_list = [fReadname, rReadname]
    outfilestring = str()
    for line in config:
        line = line.rstrip()
        line = line.replace(": ", "=")
        if line.startswith("max_rd_len"):
            outfilestring += line + "\n" + "[LIB]" + "\n"
        elif line.startswith("avg_ins"):
            outfilestring += line + "\n"
        elif line.startswith("reverse_seq"):
            outfilestring += line + "\n"
        elif line.startswith("asm_flags"):
            outfilestring += line + "\n"
        elif line.startswith("rd_len_cutoff"):
            outfilestring += line + "\n"
        elif line.startswith("rank"):
            outfilestring += line + "\n"
        elif line.startswith("pair_num_cutoff"):
            outfilestring += line + "\n"
        elif line.startswith("map_len"):
            outfilestring += line + "\n"
    outfilestring += "q1={}\nq2={}".format(read_list[0], read_list[1])
    outfile.write(outfilestring)
