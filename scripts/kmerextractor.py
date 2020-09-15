import sys

kmerlogfile = sys.argv[1]
lower_file = sys.argv[2]
upper_file = sys.argv[3]
optimal_file = sys.argv[4]

with open(kmerlogfile, "r") as kmerlog:
    for line in kmerlog:
        line = line.rstrip()
        if line.startswith("best k:"):
            best_k = int(line.split(": ")[1])
            lower_k = best_k - 4 
            upper_k = best_k + 4
with open(lower_file, "w") as lower:
    lower.write(str(lower_k))
with open(upper_file, "w") as upper:
    upper.write(str(upper_k))
with open(optimal_file, "w") as optimal:
    optimal.write(str(best_k))
