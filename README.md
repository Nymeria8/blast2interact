# blast2cytoscape

This program takes the Blastx output in tabular form against a database with interactions between sequences known (file with pairs of interactions and type of interactions) and retrive
the interactions between our sequences found by homologie. 

Also recives DESEQ2 output files and organize it in other file with the names and respective fold changes

This two files can be directly feeded to cytoscape to network construction (http://www.cytoscape.org/)


##Usage:



In the command line:
<pre><code>python3 python interactions.py blastoutput correspondance_file output_correspondance output_expression DESEQ2_ouput_files(variable)
</code></pre>


Notes:

+ This script was developed to use with the file protein.actions and respective sequences from the string database (http://string-db.org/) 
but it can be used with other databases or files of interactions, since the file is formated as gene_one {tab} gene_two {tab} interaction, and the names in the files
corresponds to the ones of the databases.
+ This script takes a blastx output, so its necessary to run it against the database of sequences with outfmt=7, or another tabular type file
+ It only takes he best hit of blast
+ The script sets a cutoff of identity > 90% and alignment>100aa, that can be change directly in the script
+ The ouput with the interactions can be write with the databases names or our sequences names (default). The script can be easily change to use one or another

##Dependencies:

+ Python3



##License:

GPLv2


##Found a bug?

Or maybe just wan to drop some feedback? Just open an issue on github!
   
