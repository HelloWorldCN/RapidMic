RapidMic is a simple, easy-to-use, rapid for computing Maximal Information-based Nonparametric Exploration(D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, E. Lander, M. Mitzenmacher, and P. Sabeti.  Detecting
novel associations in large datasets. Science, 6062(334):1518¨C1524, 2011.)
In this implementation , we use optimized multiple threads to speeded up processing based on Davide Albanese et al.'s code at http://minepy.sourceforge.net/.
This software can help to identify interesting relationships between pairs of variables in large data sets is increasingly important.
This document explains the use of RapidMic.

RapidMic is available at https://github.com/HelloWorldCN/RapidMic



#Table of Contents
=================

* Quick Start
* Installation and Data Format


#Quick Start
============================
### Usage: 
`RapidMic -i datasetfile [inputoptions] [outputoptions]`

### options:
	-h, --help     display this help and exit;
	-V, --version     output version information and exit;
	-i <file>, --input=<file>     input filename;
	-l <label style>, --inputlabel=<label style>     if input csv file's first line is comlumn name and the first comlumn is row name, then <label style>=0 ,else if just only comlumn name <label style>=1, else if just only row name <label style>=2 ;
	-a <alpha value>, --alpha=<alpha>     the exponent in B(n) = n^alpha (default: 0.6.) alpha must be in (0, 1.0]
	-c <clumps value>, --clumps=<c>     determines how many more clumps there will be than columns in every partition. Default value is 15,c must be > 0;
	-o <file>, --output=<file>     output filename (default: mine_out.csv);
	-L <label style>, --outputlabel=<label style>     output csv file adopt number index of row as vaiable label ,then <label style>=0£¬else adopt input file's row name, then <label style>=1;
	-A <allPairs>, --allPairs     will cause MINE to compare all pairs of variables against each other; 
	-b <var index>, --pairsBetween=<var index>     will compare each of the first i variables to each of the rest of the variables. Variables are indexed from 0;input variable <var index> must be in (0, number of variables in file)
	-m <var index>, --master=<var index>     variable <var index> vs. all <var index> must be in [0, number of variables in file);
	-p <var1 index> <var2 index>,--onePair=<var1 index> <var2 index>     variable <var1 index> vs variable <var2 index> <var1 index> and <var2 index> must be in [0, number of variables in file);	



If output file is omitted, the compute result will be printed on the screen.


Example
=======
First download datasets from http://www.exploredata.net/Downloads

1)Ex1: Compute one pair(1 vs. 2) variables

	`RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -p 1,2 -l 0 -L 1`

2)Ex2: Compute all pairs of variables against each other

	`RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -A -l 0 -L 1`

3)Ex3:Compute the 4th variable vs. all the other

	`RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -m 4 -l 0 -L 1`

4)Ex4:Compute each of the first 4 variables to each of the rest of the variables

	`RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -b 4 -l 0 -L 1`



#Installation and Data Format
============================
You can use pre compiled binary in sub directory bin,
Or compile source yourself.
On Unix ,linux and macos systems, type the following to build in terminal 
	`g++ main.cpp core.c mine.c cppmine.cpp arg_parser.cc stringenc.cpp -o RapidMic -lpthread`

On Windows, you can use VS2010 to build

### input file format
The format of input data file is comma-separated values (CSV) file:
### 
		name,colname1,colname2,...
		variable1,1,2,...
		variable2,3,5,...
		variable3,6,9,...
		variable4,...
		...
		

Each line contains an variable instance and is ended by a '\n' character. 


### output file format
The format of output data file is comma-separated values (CSV) file:
###
		var1,var2,mic,mev,mcn,mas
		0,1,0.1,0.2,0.4,0.5
		....
	
###or
		var1,var2,mic,mev,mcn,mas
		YAL001C,YAL014C,0.1,0.2,0.4,0.5
		....
	
Each line contains one pair's compute result and is ended by a '\n' character. 	
	

 