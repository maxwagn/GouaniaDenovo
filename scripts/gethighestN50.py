import sys
import pandas as pd


####### IMPORTANT CHANGE LATER 
####### "soapoutfile.write('../assemblies/{}_SOAPDENOVO/{}_soap_{}.scafSeq'.format(ID, ID, SOAPmaxN50))"
####### TO
####### "soapoutfile.write('assemblies/{}_SOAPDENOVO/{}_soap_{}.scafSeq'.format(ID, ID, SOAPmaxN50))"

# usage: python gethighestN50.py <ID> <metadata> <quast_report_transposed> <outputdir>
ID = sys.argv[1]
metadata = sys.argv[2]
quast = sys.argv[3] ### use transposed table
outputdir = sys.argv[4]

mt = pd.read_csv(metadata, dtype=str, sep='\t')
id_list = mt["sequence_id"].tolist()
report = pd.read_csv(quast, sep='\t')
software = ["soap", "abyss", "discovar"]
ksize = ["lowerK", "optimalK", "upperK"]

## Prepare dataframe - add collumns software, ID and ksize for subsequent subsetting
for i in software:
    report.loc[report['Assembly'].str.contains(i),"software"] = i
for i in id_list:
    report.loc[report['Assembly'].str.contains(i),"ID"] = i
for i in ksize:
    report.loc[report['Assembly'].str.contains(i),"ksize"] = i

## Subset data and write to different files

with open('{}{}_SOAP_highestN50.txt'.format(outputdir, ID), "w") as soapoutfile:
    SOAPDataFrame = report[report['software'] == 'soap']
    SOAPDataFrame = SOAPDataFrame[SOAPDataFrame["ID"] == ID]
    SOAPDataFrame = SOAPDataFrame.set_index("ksize")
    SOAPmaxN50 = SOAPDataFrame['N50'].idxmax()
    soapoutfile.write('../assemblies/{}_SOAPDENOVO/{}_soap_{}.scafSeq'.format(ID, ID, SOAPmaxN50))

with open('{}{}_ABYSS_highestN50.txt'.format(outputdir, ID), "w") as abyssoutfile:
    AbyssDataFrame = report[report['software'] == 'abyss']
    AbyssDataFrame = AbyssDataFrame[AbyssDataFrame["ID"] == ID]
    AbyssDataFrame = AbyssDataFrame.set_index("ksize")
    AbyssmaxN50 = AbyssDataFrame['N50'].idxmax()
    abyssoutfile.write("../assemblies/{}_ABySS/{}_abyss_{}-scaffolds.fa".format(ID, ID, AbyssmaxN50))

with open('{}{}_DISCOVAR_highestN50.txt'.format(outputdir, ID), "w") as discovaroutfile:
    DiscovarDataFrame = report[report['software'] == 'discovar']
    discovaroutfile.write("../assemblies/{}_DISCOVAR/a.final/{}_discovar.fasta".format(ID, ID))
