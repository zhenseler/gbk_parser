# gbk_parser
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
