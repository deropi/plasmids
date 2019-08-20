import os
import itertools
import sys
import argparse
from Bio import SeqIO
##### Argument parsing
parser = argparse.ArgumentParser(description = "This script analyzes assembled genomes and identifies plasmis present in RefSeq", epilog = 'Please, report bugs to dpicazo@ifam.uni-kiel.de')
parser.add_argument("contigs",metavar="contigs",type=str,help="Fasta file containing contigs to parse")
parser.add_argument("coverage",metavar="coverage",type=float,help="Minimum query coverage allowed (in percentage)")
parser.add_argument("identity",metavar="identity",type=float,help="Minimum sequence identity")
parser.add_argument('output_name', metavar='output_name', type=str, help='Provide with an output name')
parser.add_argument("-db", "--database", help="Specify the path to the plasmids database. If not specified, it will download the current version and create blast databases",action="store")
args=parser.parse_args()
##############################################################
# 1. Here, we retrieve current plasmid sequences from NCBI and generate blast database
#############################################################
def download_db():
	plasmids=r'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/*genomic.fna*'
	print 'Creating plasmids database...'
	os.system("mkdir plasmids_db")
	os.system("wget "+plasmids + " -P plasmids_db")
	for db in os.listdir('plasmids_db'):
	        os.system('gzip -d ./plasmids_db/' + db)
        	db = './plasmids_db/' +  db[:-3]
	        os.system('makeblastdb -in %s -dbtype nucl' %db)


if not args.database:
	download_db()
	plasmids_path = 'plasmids_db'
else:
	plasmids_path = args.database

if not plasmids_path.endswith('/'):
	plasmids_path += '/'
#########################################################
# 2. Now, let's store the contigs in a dictionary
#######################################################

def read_fasta(file, dict):
	header = ''
	with open(file, 'r') as fasta:
		for line in fasta:
			if line.startswith('>'):
				if header != '':
					dict[header] = seq
				header = line.strip().split()[0][1:]
				seq = ''
			else:
				seq += line.strip()
	dict[header] = seq

contigs = {}
read_fasta(args.contigs, contigs)

queries = {}

#print contigs.keys()
## that was totally useless
## Here, we blast the contigs against plasmids database
for db in [file for file in os.listdir(plasmids_path) if file.endswith('.fna')]:
	os.system('blastn -db ' + plasmids_path + db + ' -query ' + args.contigs + ' -outfmt 7 -evalue 10e-10 -num_threads 10 -out ' + plasmids_path + db[:-3]+ args.output_name + '.blastn')
	read_fasta(plasmids_path + db, queries)
	
## Let's make a little filtering...
hits = {}
plasmids_hits = {}
plasmids_hit_lengths = {}

blast_clean = open(args.output_name+'_blastn.clean', 'w')
plasmids_file = open(args.output_name+'_plasmids_sum.tab', 'w')

for blastn in [file for file in os.listdir(plasmids_path) if file.endswith(args.output_name+'.blastn')]:
	with open(plasmids_path+ blastn, 'r') as input:
		for line in input:
			if not line.startswith('#'):
				info = line.strip().split()
				if not hits.has_key(info[0]):
					query  = info[0]
					subject = info[1]
					id = float(info[2])
					query_cove = 100*(float(info[3])/float(len(contigs[query])))
					if id >= args.identity and query_cove >= args.coverage:
						hits[query] = info[1:]
						blast_clean.write(line)
						if not plasmids_hits.has_key(subject):
							plasmids_hits[subject] = []
							plasmids_hit_lengths[subject] =set()
							#plasmids[subject][query] = set()
						plasmids_hits[subject].append(query)
						if int(info[8]) < int(info[9]):
							for i in range(int(info[8]), int(info[9])):
								plasmids_hit_lengths[subject].add(i)
						else:
							for i in range(int(info[9]), int(info[8])):
                                                                plasmids_hit_lengths[subject].add(i)
#print plasmids
blast_clean.close()

plasmids_file.write('Plasmid_ID\tplasmid_coverage\tmatching_contigs\tlength_plasmid\n')
for plasmid in plasmids_hits:
#	print plasmids[plasmid].values()
	plasmid_covered = 100*(float(len(plasmids_hit_lengths[plasmid]))/float(len(queries[plasmid])))
	plasmids_file.write(plasmid + '\t' + str(plasmid_covered) + '\t' + str(plasmids_hits[plasmid]) + '\t' + str(len(queries[plasmid]))  + '\n'  )
plasmids_file.close()


