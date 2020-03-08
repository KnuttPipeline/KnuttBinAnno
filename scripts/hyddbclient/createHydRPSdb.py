#!/usr/bin/env/python3

import os
import sys
import csv
import urllib.request
import argparse
import subprocess

superfamilyfile_url = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links"
pssms_url = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz"
idlookup_url = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz"
superfamilies = {"NiFe":"cl21493"}
pssms = {"Fe":{"pfam03201"},"NiFe":set(),"FeFe":{"pfam02906","COG4624","TIGR02512"}}


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser(description="Download the hydrogenase PSSMs from NCBI CDS and create a RPS database to use it for finding Hydrogenase candidates. Afterwards you can run rpsblast like this:\nrpsblast -i CDS.faa -p T -e 0.00001 -m 8 -o hyd_rpsblast.tsv -d OUTDIR/DB_NAME")
parser.add_argument("outdir", type=dir_path,default=".",help="The directory to store the downloaded models")
parser.add_argument("-d","--db_name",default="hyd",help="Name of the database to create, will be stored in the outdir. Default is hyd")

args = parser.parse_args()

targetdir = args.outdir
dbname = args.db_name

idlookup_file = os.path.join(targetdir,"cddid_all.tbl.gz")
if not os.path.isfile(idlookup_file):
	print("Downloading id lookup table")
	urllib.request.urlretrieve(idlookup_url,idlookup_file)



superfamilyfile = os.path.join(targetdir,"family_superfamily_links.tsv")
superfamily_index = 2
pssm_index = 0

if not os.path.isfile(superfamilyfile):
	print("Downloading superfamily table")
	urllib.request.urlretrieve(superfamilyfile_url, superfamilyfile)


with open(superfamilyfile,"r") as superfamilies_stream:
	reader = csv.reader(superfamilies_stream, delimiter='\t')
	for line in reader:
		if len(line) >= 2:
			for center,superfamily in superfamilies.items():
				if line[superfamily_index] == superfamily:
					pssms[center].add(line[pssm_index])
print(pssms)


pssm_files = [pssm+".smp" for pssmset in pssms.values() for pssm in pssmset ]
print(pssm_files)

pssms_file = os.path.join(targetdir,"cdd.tar.gz")

if not os.path.isfile(pssms_file):
	print("Downloading models ...")
	urllib.request.urlretrieve(pssms_url,pssms_file)


print("Extracting selected models ...")
exit_code = subprocess.call(["tar","-xzvf",pssms_file,"-C",targetdir]+pssm_files)
if exit_code != 0:
	print ("Error extracting files!")
	sys.exit(1)


pnfile_path = dbname+".pn"
with open(os.path.join(targetdir,pnfile_path),"w") as pnfile:
	for pssm_file in pssm_files:
		pnfile.write(pssm_file)
		pnfile.write("\n")

tsvfile_path = dbname+".tsv"
with open(os.path.join(targetdir,tsvfile_path),"w") as tsvfile:
	tsvfile.write("center\tpssm")
	tsvfile.write("\n")
	for center,pssms in pssms.items():
		for pssm in pssms:
			tsvfile.write(f"{center}\t{pssm}")
			tsvfile.write("\n")


print("Creating database...")
os.chdir(targetdir)
exit_code = subprocess.call(["makeprofiledb","-in",pnfile_path,"-out",dbname])
if exit_code != 0:
	print ("Error building db!!")
	sys.exit(1)
print("Done!")
