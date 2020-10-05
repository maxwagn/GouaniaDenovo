import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

BUSCOALLINALL = sys.argv[1]
quast = sys.argv[2] #quast report transposed_report.tsv
metadata = sys.argv[3]

mt = pd.read_csv(metadata, dtype=str, sep='\t')
id_list = mt["sequence_id"].tolist()
software = ["SOAP", "ABYSS", "DISCOVAR"]


buscostats = ["C","S","D","F","M","n"]

BUSCO_csv = str()

with open(BUSCOALLINALL, "r") as busco:
    for line in busco:
        line = line.rstrip()
        if line.startswith("# Summarized"):
            seqID = line.split("/")[-1]
            seqID = seqID.split(".fasta")[0] 
        if line.startswith("\tC:"):
            stats = line.replace("[",",")
            stats = stats.replace("]","")
            stats = stats.replace("\t","")
            stats = stats.replace("%","")
            for i in buscostats:
                stats = stats.replace("{}:".format(i),"")
            BUSCO_csv += "{},{}\n".format(seqID, stats)

BUSCO_csv = BUSCO_csv[:-1]

busco_report = pd.DataFrame([x.split(',') for x in BUSCO_csv.split('\n')])
busco_report.columns = ["Assembly","Complete","Single","Duplicated","Fragmented","Missing","Total"]

## Prepare dataframe - add collumns software, ID and ksize for subsequent subsetting
for i in software:
    busco_report.loc[busco_report['Assembly'].str.contains(i),"software"] = i
for i in id_list:
    busco_report.loc[busco_report['Assembly'].str.contains(i),"sequence_id"] = i

quast_report = pd.read_csv(quast,sep='\t',header=(0))

final_report_raw = pd.merge(busco_report, quast_report, on = ["Assembly"])
final_report_raw = final_report_raw.rename(columns={"# N's per 100 kbp": "NsPer100kbp"})
final_report_raw.to_csv("../reports/final/BUSCOandQUAST_summary_final.tsv", sep="\t")

graphsIDs = ["{}_{}".format(i,j) for i, j in zip(final_report_raw["software"], final_report_raw["sequence_id"])] 
sing = final_report_raw["Single"].to_numpy(dtype='float')
dup = final_report_raw["Duplicated"].to_numpy(dtype='float')
frag = final_report_raw["Fragmented"].to_numpy(dtype='float')
miss = final_report_raw["Missing"].to_numpy(dtype='float')
width = 0.4


busco_plot = plt.figure(figsize=(8,6))
plt.bar(graphsIDs, sing, label='Single')
plt.bar(graphsIDs, dup, bottom=sing, label='Duplicated')
plt.bar(graphsIDs, frag, bottom=dup+sing, label='Fragmented')
plt.bar(graphsIDs, miss, bottom=dup+sing+frag, label='Fragmented')
plt.xticks(rotation=90)
plt.title('% BUSCO genes', size = 16)
plt.ylabel('% BUSCO genes', size = 14)
plt.legend(loc='center', bbox_to_anchor=(1.15, 0.5))
plt.tight_layout()
plt.savefig('../reports/final/BUSCO_summary_PLOT.pdf')

quastPlot = plt.figure(figsize=(12,6)) # Create matplotlib figure
ax = quastPlot.add_subplot(111) # Create matplotlib axes
ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax.
plt.title("Summary - N50 & N's per 100 kbp", size = 16)
final_report_raw.NsPer100kbp.plot(kind='bar', color="C9", ax=ax2, width=width, position=0)
final_report_raw.N50.plot(kind='bar', color="orange", ax=ax, width=width, position=1)
ax.set_xticklabels(graphsIDs)
ax.set_ylabel("N50", color = "orange", size = 14)
ax.tick_params(axis='y', colors='orange')
ax2.set_ylabel("# N's per 100 kbp", color = "C9", size = 14)
ax2.tick_params(axis='y', colors='C9')
plt.xlim((-0.8, len(final_report_raw)-0.2)) # makes edges of graph
plt.tight_layout()
plt.savefig('../reports/final/QUAST_summary_PLOT.pdf')
