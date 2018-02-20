'''
Obtain DNA sequences for features containing input tag in .gbk file


Dependencies:

	- None
	
Usage:

	- python gbk_parser.py <input_gbk> '<tag>' <output_folder>
	
		- input_gbk:		Your .gbk file
		
		- tag:				Script will get sequences from features containing
							the input tag in their name
						
		- output_folder:	Folder resulting fasta file will be written to
		
	- Have a bunch of .gbk files and want to search for same tag in all?
	
			in Unix shell...
			
		- for file in <path/to/files>*.gbk; do python gbk_parser.py $file '<tag>' <output_folder>; done
		
			and if you want to combine all the files in the end...
			
		- cat <output_folder>/*_<tag>_parsed.gbk > <combined_filename>

	
Workflow:

	1. Grab coords from lines with '..'
	2. Add coords to list if they contain a feature containing input tag
	3. Grab DNA sequence for each coord pair in coords list
	4. Write out all unique instances of feature
	
To-do:
	
	- Work in ability to blast sequences
	- Automatically decide to perform on one file, or if pointed to a folder, all
	  files in the folder and cat them together at the end
	- Test in 3.X

'''



# imports

from sys import argv
import re
from subprocess import call			# Unix calls in Python


# args

file_name = argv[1] 			# Input file name
goi = str(argv[2]) 				# Must be a string identical to what you
								# would find in the annotation file
output_folder = argv[3]			# Folder you want output in


# File names and such

output_prefix = file_name[:-4]
sample_name = output_prefix.split('/')[-1]

fasta_prefix = '>'+sample_name

goi_output = goi.replace(' ','_')

seq_file = open(goi_output+'_seq_file.fasta','a')

call('mkdir '+output_folder,shell=True)
if output_folder[-1] != '/':
	output_folder += '/'
output_str = '%s%s_%s_parsed.fasta' % (output_folder,sample_name,goi_output)
fasta_output = open(output_str,'w')


# Initialize before reference

coord_line = 0
coord_list = []
contig = ''
seq_list = []


# Grab line, coordinates for each goi instance

for i,l in enumerate(open(file_name,'U')):

	if '..' in l:
		coord_line = i
		type = l.split('     ')[1]
		current_coords = re.findall('\d+',l)
		if 'complement' in l:
			comp = True
		else:
			comp = False
		start_coord = current_coords[0]
		end_coord = current_coords[1]

	elif goi in l:
		fasta_suffix = l.split('"')[1]				# Grab goi substring
		coord_tup = ((fasta_suffix,int(coord_line),int(start_coord),int(end_coord),comp))
		coord_list.append(coord_tup)


# Grab DNA sequence for each coordinate

for coord_tup in coord_list:
	for i,l in enumerate(open(file_name,'U')):
		if i > coord_tup[1]:
			if l[:2] == '//':
				break
			elif l.rstrip() == 'ORIGIN':
				origin_line = i
			else:
				if l.strip()[0] in str([1,2,3,4,5,6,7,8,9]):
					contig = contig + l.replace(' ','')
	contig = ''.join([i for i in contig if not i.isdigit()])
	contig = contig.replace('\n','')
	seq = contig[coord_tup[2]-1:coord_tup[3]]			# Slice contig for goi sequence
	if comp == True:
		seq = rev_comp(seq)								# rev_comp if needed
	seq_tup = (seq,origin_line)
	if not seq_tup in seq_list:							# Prevent adding same sequence multiple times
		seq_list.append(seq_tup)
		
for i,tup in enumerate(seq_list):
	fasta_output.write('%s_%s_%s\n%s\n' % (sample_name,coord_tup[0].replace(' ','_'),str(i+1),str(seq)))


# If goi not in .gbk file

if seq_list == []:
	print('\nNo instance of %s in %s.gbk\n' % (goi,sample_name))


fasta_output.close()



# Function from bioinfo_tools

def rev_comp(input_seq):
	rev_comp_dict = {'A':'T','C':'G','G':'C','T':'A'}
	return ''.join(rev_comp_dict[nuc] for nuc in input_seq[::-1])







