#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from requests_toolbelt.multipart.encoder import MultipartEncoder
from bs4 import BeautifulSoup


import os
import argparse
import sys
import requests
import re
import json
import shutil
from shlex import quote
import subprocess
from time import sleep
import tempfile
import io
import operator
import functools
import itertools


def dir_path(parser,arg):
    if not os.path.isdir(arg):
        parser.error(f"The directory {arg} does not exist!")
    else:
       return arg

def file_path(parser,arg):
    if not os.path.isfile(arg):
       parser.error(f"The file {arg} does not exist!")
    else:
       return arg


parser = argparse.ArgumentParser(description="Run RPSblast on the given sequences to find hydrogenases and  classify those with the webservice Hyddb. Remeber to cite them in your work.")
parser.add_argument("input", type=lambda x: file_path(parser, x),help="FASTA formatted protein sequences. Default is to read from STDIN.")
parser.add_argument("-o","--output", nargs="?", type=argparse.FileType("w"), default=sys.stdout,help="Location of the result file. Default is STDOUT")
parser.add_argument("-l","--log", nargs="?", type=argparse.FileType("w"), default=sys.stderr,help="location of the log output. Default is STDERR")
parser.add_argument("-f","--dbdir", type=lambda x: dir_path(parser, x),default=".",help="The directory where the database was created with the createHydRPSdb.py script.")
parser.add_argument("-n","--db_name",default="hyd",help="Name of the database. Default is hyd")
parser.add_argument("-d","--delay",nargs="?",type=int,default=10,help="number of seconds to wait between each status request")
parser.add_argument("-r","--retries",nargs="?",type=int,default=200,help="number of times to ask for status updates before just failing, 0 means unlimited retries. Default is ")
parser.add_argument("--chunksize",nargs="?",type=int,default=25,help="The number of sequences to classify in one go. Default is 25")
parser.add_argument("-e","--rpseval",type=float,default=0.00001,help="E-value for the RPS Blast. Default is 10^-5")



##
parser.add_argument("--contigregex",type=str,default=r".*cid=([^ ]+)(?: |$)",help="The regex of the field with contig id (>cds_id foo=bar cid=k127_7471)")
##


baseurl = "https://services.birc.au.dk/hyddb/"


def parseRPSdata(rps,center_file,asc_file):
	hitdata = pd.read_csv(rps,names=["query","subject","identity","allength","mismatches","gapopenings","querystart","queryend","subjectstart","subjectend","eval","score"],sep="\t")
	#print(hitdata)
	if hitdata.empty:
		return None
	hitdata["asc"] = hitdata.subject.str.extract(r"(\d+)$").astype("int64")
	hitdata = hitdata.set_index("asc")

	centerdata = pd.read_csv(center_file,sep="\t",index_col="pssm")
	ascdata = pd.read_csv(asc_file,names=["asc","pssm","shortname","name","subjectlength"],compression="gzip",sep="\t",index_col="pssm")
	centerdata = centerdata.join(ascdata).set_index("asc")

	data = hitdata.join(centerdata)
	data = data.set_index("query")
	return data

def callRPS(proteinfastapath,rpseval,dbdir,dbname,logstream):
	db = os.path.join(dbdir,dbname)
	center_file = os.path.join(dbdir,dbname+".tsv")
	asc_file = os.path.join(dbdir,"cddid_all.tbl.gz")
	with tempfile.NamedTemporaryFile() as tmp: #-query -evalue -outfmt -db
		cmd = ['rpsblast', '-query',quote(proteinfastapath),"-evalue",str(rpseval),"-outfmt","6","-db",quote(db)]
		print(" ".join(cmd),file=logstream)
		returncode = subprocess.call(cmd,stdout=tmp, stderr=logstream)
		if returncode != 0:
			tmp.seek(0)
			#print(tmp.read().decode(sys.stdout.encoding),file=logstream)
			raise subprocess.CalledProcessError(returncode,cmd[0])
		tmp.seek(0)
		data = parseRPSdata(tmp,center_file,asc_file)
		return data
	return None


def selectLongestNonOverlappingEntries(dataframe):
	dataframe = dataframe.reset_index().sort_values(by="querystart",ascending=True)
	selectedrows = [0]
	doesNotoverlap = lambda row1,row2: not (row1["queryend"] >= row2["querystart"] and row2["queryend"] >= row1["querystart"])
	if(len(dataframe.index)>1):
		for rowi in range(1,len(dataframe.index)):
			if all([doesNotoverlap(dataframe.iloc[rowi],dataframe.iloc[lasti]) for lasti in selectedrows ]):
				selectedrows.append(rowi)
	return dataframe.iloc[selectedrows]

def getCRSFToken(session,url=baseurl):
  classifypage = session.get(url)
  classifypage.raise_for_status()
  soup = BeautifulSoup(classifypage.content,features="html.parser")
  tag = soup.find("input",attrs={"name":"csrfmiddlewaretoken"})
  # TODO fail on missing tag or attribute 
  return tag["value"]

def getJobId(content):
	soup = BeautifulSoup(content,features="html.parser")
	link = soup.find(string="Link to this page").findNext("a").attrs["href"]
	jobid = re.search(r"(?<=/jobs/).+(?=/wait)",link).group(0)
	return jobid

def submitJob(session,instream):
  token = getCRSFToken(session)
  m = MultipartEncoder(fields={"csrfmiddlewaretoken":token,"sequences_file":("sequences.fasta",instream,"text/plain")})
  result = session.post(baseurl,data=m,headers={'Content-Type': m.content_type})
  result.raise_for_status()
  # TODO error checking
  jobid = getJobId(result.content)
  return jobid

def getJobStatus(session,jobid):
  statusurl=f"{baseurl}jobs/{jobid}/status"
  result = session.get(statusurl)
  result.raise_for_status()
  status = json.loads(result.content)["status"]
  return status

def downloadResults(session,jobid):
  resulturl=f"{baseurl}jobs/{jobid}/results.csv"
  data = pd.read_csv(resulturl,sep=";",names=["hitid","Classification"],index_col="hitid")
  return data

def reClassifyWithDownStream(session,jobid,hydwithdownstream_stream):
	resultpageurl=f"{baseurl}jobs/{jobid}/results"
	resultpage = session.get(resultpageurl)
	resultpage.raise_for_status()
	soup = BeautifulSoup(resultpage.content,features="html.parser")
	additionalre = re.compile(r"^id_custom_(.+)$")
	textareas = soup.find_all("textarea",id=additionalre)
	if len(textareas) == 0:
		return None
	sequencefields = {additionalre.match(textarea.attrs["id"]).group(1).strip():textarea.attrs["name"] for textarea in textareas}
	token = textareas[0].find_parent("form").findChild("input",attrs={"name":"csrfmiddlewaretoken"}).attrs["value"]
	fieldstosend = {"csrfmiddlewaretoken":token}
	for record in SeqIO.parse(hydwithdownstream_stream, "fasta"):
		nameforseq = sequencefields.get(record.id)
		if nameforseq:
			fieldstosend[nameforseq]=str(record.seq)
	newclass_result = session.post(resultpageurl, data = fieldstosend)
	newclass_result.raise_for_status()
	jobid = getJobId(newclass_result.content)
	return jobid

def waitForJob(session,jobid,lprint):
	status = getJobStatus(session,jobid)
	currenttry = 0
	while currenttry<args.retries and status != "SUCCESS" and status!="FAILURE":
		lprint(f"Last Status: {status}")
		sleep(args.delay)
		status = getJobStatus(session,jobid)
		currenttry+=1
	return status

def main(args,logfile,outfile):
	lprint = lambda *pargs:print(*pargs,file=logfile, flush=True)
	lprint("Started Hyddb Client...")
	lprint("Arguments: "+str(args))
	session = requests.Session()
	lprint("Running RPS ...")
	data = callRPS(args.input,args.rpseval,args.dbdir,args.db_name,logfile)
	#print(data)
	if data is None or data.empty:
		_ = outfile.write("query	subject	identity	allength	mismatches	gapopenings	querystart	queryend	subjectstart	subjectend	eval	score	center	shortname	name	subjectlength	hitid	Classification\n")
		return
	lprint(f"{len(data.index)}  Hit(s)!")
	lprint("Extracting longest non overlapping hits for every query and center ...")
	entries = data.groupby(["query","center"]).apply(selectLongestNonOverlappingEntries).reset_index(drop=True)
	lprint(f"{len(entries.index)}  Hit(s) remaining!")
	lprint("Extracting hit regions from query file ...")

	def contiggetter(record): return re.match(args.contigregex,record.description).group(1)
	records = list(SeqIO.parse(args.input, "fasta"))
	records.sort(key=contiggetter)
	records={contig:list(sorted(g,key=lambda record: record.id)) for contig,g in itertools.groupby(records, key=contiggetter)}

	sequences = []
	for contig,recordgroup in records.items():
		for i,record in enumerate(recordgroup):
			entries_for_record = entries[entries["query"]==record.id]
			if len(entries_for_record.index)==0:
				continue
			for index,row in entries_for_record.iterrows():
				newid = f"{record.id}:{row.querystart}:{row.queryend}"
				hyd_andds_sequence = SeqIO.SeqRecord(record.seq[row.querystart:],id=newid,description="")
				if i!=len(recordgroup)-1:
					hyd_andds_sequence+=functools.reduce(operator.add,recordgroup[i+1:])
				hyd_andds_sequence.id = newid
				hyd_andds_sequence.description = ""
				hyd_sequence = SeqIO.SeqRecord(record.seq[row.querystart:row.queryend],id=newid,description="")
				sequences.append((hyd_sequence,hyd_andds_sequence))
	chunksize = args.chunksize
	sequences = [sequences[i * chunksize:(i + 1) * chunksize] for i in range((len(sequences) + chunksize - 1) // chunksize )]
	lprint(f"Hits split into {len(sequences)} chunk(s)!")
	curchunk = 1
	results = []
	#print(list(map(len,sequences.values())))
	for sequences_chunk in sequences:
		with tempfile.NamedTemporaryFile() as tmp:
			lprint(f"Saving chunk {curchunk}/{len(sequences)} into temp file ...")
			tmp_wrapped = io.TextIOWrapper(tmp,"UTF-8")
			SeqIO.write([sequence_pair[0] for sequence_pair in sequences_chunk],tmp_wrapped,"fasta")
			lprint([sequence_pair[0].id for sequence_pair in sequences_chunk])
			tmp_wrapped.flush()
			lprint("Submitting job...")
			tmp.seek(0)
			jobid = submitJob(session,tmp)
			lprint(f"Got jobid: {jobid}")
			status = waitForJob(session,jobid,lprint)
			if status=="SUCCESS":
				lprint("Retrieving results...")
				newresults = downloadResults(session,jobid)
				lprint("Saving hydrogenases with their downstream sequences for further classification ...")
				tmp_wrapped.seek(0)
				SeqIO.write([sequence_pair[1] for sequence_pair in sequences_chunk],tmp_wrapped,"fasta")
				tmp_wrapped.flush()
				tmp_wrapped.seek(0)
				additionaljobid = reClassifyWithDownStream(session,jobid,tmp_wrapped)
				if additionaljobid:
					lprint(f"Got further jobid: {additionaljobid}")
					lprint("Waiting for further classification ...")
					status = waitForJob(session,additionaljobid,lprint)
					if status=="SUCCESS":
						lprint("Retrieving further results...")
						additionalresults = downloadResults(session,additionaljobid)
						newresults.update(additionalresults)
					else:
						lprint("Failed!")
						sys.exit(1)
				lprint(newresults)
				results.append(newresults.reset_index(drop=False))
			else:
				lprint("Failed!")
				sys.exit(1)
		curchunk = curchunk + 1
	lprint("Merging results...")
	results = pd.concat(results,ignore_index=True)
	results["hitid"] = results["hitid"].apply(lambda val: val.strip())
	entries["hitid"] = entries.apply(lambda row: f"{row.query}:{row.querystart}:{row.queryend}",axis=1)
	results = entries.merge(results,"left",on="hitid")
	
	results.to_csv(outfile,sep="\t",header=True,index=False)


if __name__== "__main__":
	args = parser.parse_args()
	main(args,args.log, args.output)
