RapidMic is a simple, easy-to-use, rapid for computing Maximal Information-based Nonparametric Exploration(D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, E. Lander, M. Mitzenmacher, and P. Sabeti.  Detecting
novel associations in large datasets. Science, 6062(334):1518-1524, 2011.)
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
	-i <file>, --input=<file>     input comma-separated value(CSV) file filename;The filename input is a string enclosed in single quotes;The file can only contain numeric values;
	-l <label style>, --inputlabel=<label style>     if input csv file's first line is comlumn name and the first comlumn is row name, then <label style>=0 ,else if just only comlumn name <label style>=1, else if just only row name <label style>=2 ;
	-a <alpha value>, --alpha=<alpha>     the exponent in B(n) = n^alpha (default: 0.6.) alpha must be in (0, 1.0]
	-c <clumps value>, --clumps=<c>     determines how many more clumps there will be than columns in every partition. Default value is 15,c must be > 0;
	-o <file>, --output=<file>     output filename (if not be set, will print result in screen);
	-L <label style>, --outputlabel=<label style>     output csv file adopt number index of row as vaiable label ,then <label style>=0 ;else adopt input file's row name, then <label style>=1;
	-A <allPairs>, --allPairs     compare all pairs of variables against each other; 
	-b <var index>, --pairsBetween=<var index>      compare each of the first i variables to each of the rest of the variables. Variables are indexed from 0;input variable <var index> must be in (0, number of variables in file)
	-m <var index>, --master=<var index>     compare variable <var index> to the rest of the variables, <var index> must be in [0, number of variables);
	-p '<var1 index>,<var2 index>',--onePair='<var1 index>,<var2 index>'    compare one pair variables <var1 index> and <var2 index>, <var1 index> and <var2 index> must be in [0, number of variables);	


### Notice:
	If output file is omitted, the compute result will be printed on the screen.
	In C++/C runtime environment, the row and column arguments are zero based, so that row=0 and col=0 specify the first .


Example
=======
First download datasets from http://www.exploredata.net/Downloads

1)Ex1: Compute one pair(1 vs. 2) variables

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -p 1,2 -l 0 -L 1

2)Ex2: Compute all pairs of variables against each other

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -A -l 0 -L 1

3)Ex3:Compute the 4th variable vs. all the other

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -m 4 -l 0 -L 1

4)Ex4:Compute each of the first 4 variables to each of the rest of the variables

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -b 4 -l 0 -L 1



#Installation and Data Format
============================
You can use pre compiled binary in sub directory bin,
Or compile source yourself.
On Unix ,linux systems and Macos type the following to build in terminal 

	g++ main.cpp core.c mine.c cppmine.cpp arg_parser.cc stringenc.cpp -o RapidMic -pthread -O3

On Macos, we recommend xcode to compile. 

On Windows, you can use VS2010 to build. If your OS is WIN64, please copy all files in the sub directory "lthread-win\dll\x64" to windows systems directory.
If your OS is WIN32, copy all files in "lthread-win\dll\x86" windows systems directory. We also provide a makefile named as ‘makefile.vc’, you can execute them with the command:

	nmake -f makefile.vc

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
	

 
 
 
=======================
RapidMic是《scicene》上提出的MIC（最大信息系数）的快速并行实现，实现代码参照了Davide Albanese等的minepy。
在实现过程中对计算过程进行了大量的并行优化，因此软件能够快速的进行大量的数据的分析计算任务，不会像MINE原始作者提供
的JAVA实现软件一样在稍微大一点的数据集进行所有序列都进行比对崩溃并且速度慢的问题。

RapidMic is available at https://github.com/HelloWorldCN/RapidMic



#内容列表
=================

* Quick Start
* Installation and Data Format


#入门
============================
### 使用方式，在控制台或者是终端中输入命令: 
`RapidMic -i datasetfile [inputoptions] [outputoptions]`

### 参数:
	-h, --help     帮助;
	-V, --version     版本;
	-i <file>, --input=<file>     输入数据文件，文件格式为CSV;
	-l <label style>, --inputlabel=<label style>     输入数据文件是否包含行名或者列名，=0则表示第一行和第一列分别为行列名，=1则表示只有列名，=2表示只有行名 ;
	-a <alpha value>, --alpha=<alpha>     the exponent in B(n) = n^alpha (default: 0.6.) alpha must be in (0, 1.0]
	-c <clumps value>, --clumps=<c>     determines how many more clumps there will be than columns in every partition. Default value is 15,c must be > 0;
	-o <file>, --output=<file>     输出数据文件，文件格式为CSV(如果没有设置则直接在屏幕上输出结果);
	-L <label style>, --outputlabel=<label style>     输出数据中采用行名还是数据的顺序索引数值来标示，则表示采用索引数值来标示 ;=1则表示采用行名;
	-A <allPairs>, --allPairs     所有数据之间两两比对; 
	-b <var index>, --pairsBetween=<var index>     前i数据和剩余数据之间两两比对，附带参数取值范围为[0, 数据个数)
	-m <var index>, --master=<var index>     指定一个数据和其它数据之间进行两两比对，附带参数取值范围为 [0, 数据个数);
	-p <var1 index> <var2 index>,--onePair=<var1 index> <var2 index>    指定两个数据之间进行比对;	



If output file is omitted, the compute result will be printed on the screen.


例子
=======
First download datasets from http://www.exploredata.net/Downloads

1)Ex1: 两个序列之间(1 vs. 2) variables

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -p 1,2 -l 0 -L 1

2)Ex2: 所有数据之间两两比对

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -A -l 0 -L 1

3)Ex3:第4个数据和其它所有数据之间进行两两比对

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -m 4 -l 0 -L 1

4)Ex4:前4个数据和后面的所有数据之间进行两两比对

	RapidMic -i "Spellman.csv" -o 2.csv -a 0.6 -c 15 -b 4 -l 0 -L 1



#安装与数据格式说明
============================
如果不想自己编译则可以直接使用bin目录下已编译的程序,
在 Unix ,linux 操作系统下采用下面的命令编译

	g++ main.cpp core.c mine.c cppmine.cpp arg_parser.cc stringenc.cpp -o RapidMic -pthread -O3
	
如果是在macos下可以用上面的命令编译，也可以用xcode打开工程文件进行编译，推荐采用xcode编译，xcode编译做的优化比直接用gcc
编译要好

在Windows下可以采用VS2010 ，2008等来编译，当然也可以在cygwin等模拟环境下编译

### 输入文件的格式
输入文件采用(CSV) 格式:
### 
		name,colname1,colname2,...
		variable1,1,2,...
		variable2,3,5,...
		variable3,6,9,...
		variable4,...
		...
		

Each line contains an variable instance and is ended by a '\n' character. 


### 输出文件的格式
输出文件采用(CSV) 格式:
###
		var1,var2,mic,mev,mcn,mas
		0,1,0.1,0.2,0.4,0.5
		....
	
###或者
		var1,var2,mic,mev,mcn,mas
		YAL001C,YAL014C,0.1,0.2,0.4,0.5
		....
	
Each line contains one pair's compute result and is ended by a '\n' character. 	
	

 
 
 
 
 
 
 

 
 
 
 
 
