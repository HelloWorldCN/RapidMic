% % Help
% Note: if your operation system is linux and unix,you may need run compile.m first
% Input Parameters
% inputcsvfile     input comma-separated value(CSV) file filename;The filename input is a string enclosed in single quotes;The file can only contain numeric values. 
% alpha            the exponent in B(n) = n^alpha (default: 0.6.) alpha must be in (0, 1.0];
% clumps           determines how many more clumps there will be than columns in every partition;
% colRowName       if input file's first line is comlumn name and the first comlumn is row name, then <label style>=0 ,else if only comlumn name then <label style>=1, else if only row name then <label style>=2 ;
% analysisStyles   Analysis task type:
%                         ='A', compare all pairs of variables against each other
%                         ='b', compare each of the first i variables to each of the rest of the variables£¬with parameter i value in [1, number of variables)
%                         ='m', compare variable i to the rest of the variables, i must be in [1, number of variables];
%                         ='p', compare one pair variables i and j, i and j must be in [1, number of variables];
% 
% params           Special parameters for each kind of analysis 

%%%%%%%%%%
% Output 
% result             The result is returned in Matrix n*6. 
% 
%       var1,var2,mic,mev,mcn,mas
% 		1,2,0.1,0.2,0.4,0.5
% 		....
%          

%%%%%%%%%%
% Output Example:
% 1	2	0.389423000000000	0.389423000000000	2.58496000000000	0.0793332000000000
% 1	3	0.296060000000000	0.296060000000000	2.58496000000000	0.0519203000000000
% 1	4	0.307326000000000	0.307326000000000	2.58496000000000	0.0731126000000000
% 1	5	0.258777000000000	0.258777000000000	2.58496000000000	0.0427501000000000
% 1	6	0.588635000000000	0.588635000000000	2.58496000000000	0.0427626000000000
% 1	7	0.316317000000000	0.316317000000000	2.58496000000000	0.0231334000000000
% 1	8	0.285196000000000	0.285196000000000	2.58496000000000	0.0409440000000000

%Notice£ºIn Matlab environment, the row and column arguments are one based, so that row=1 and col=1 specify the first value in the file. The input data file and the compiled executable RapidMic need to be in the same directory for this matlab wrapper to work properly. 


%Demo, given the file Spellman.csv that contains the comma-separated values

% EX1
% result=mine_matlab('Spellman.csv',0.6,15,0,'m',1);

% EX2£º
% result=mine_matlab('Spellman.csv',0.6,15,0,'p',[3 4]);

% EX3£º
% result=mine_matlab('Spellman.csv',0.6,15,0,'b',10);

% EX4£º
% result=mine_matlab('Spellman.csv',0.6,15,0,'A');




function [result] = mine_matlab(inputcsvfile,alpha,clumps,colRowName,analysisStyles,params)
if nargin==0, % display help
    fprintf(' code\n\n');
    
    return;
end;

% Handle arguments to function
styleparam='';
if nargin<5 error('Too few input arguments');
else
    if (colRowName>2|| colRowName<0)
        error('If inout CSV file has row or column name, you  should set argument ''colRowName''(value:0-2)');
        return;
    end;
    switch analysisStyles
        case 'A'
            if nargin==5
                styleparam=' -A ';
                disp('Method is All V.S ALL')
            else
                return;
            end;
        case 'm'
            if nargin==6 &&size(params,2)==1
                
                styleparam=[' -m ',num2str(params(1))];
                disp(['Master vairable ',num2str(params(1)),' V.S Others']);
            else
                return;
                
            end
        case 'b'
            if nargin==6 &&size(params,2)==1
                
                styleparam=[' -b ',num2str(params(1))];
                disp(['First ', num2str(params(1)),' variables V.S rest of the variables'])
            else
                return;
                
            end;
        case 'p'
            if nargin==6 &&size(params,2)==2
                
                styleparam=[' -p ',num2str(params(1)),',',num2str(params(2))];
                disp(['Variable ',num2str(params(1)),' V.S ',num2str(params(2))])
            else
                return;
            end;
        otherwise
            disp('Unknown method.')
    end
    binpath='';
    if isunix
        z='/';
         binpath=[pwd,z,'RapidMic'];
    else
        z='\';
    end;
    
    if ismac
        binpath=[pwd,z,'RapidMic-macos'];
    end;
    if ispc
        binpath=[pwd,z,'RapidMic.exe'];
    end;
        
    command =[ binpath,' -i "',inputcsvfile,'" -o "temp_out.csv"',' -a ',num2str(alpha),' -c ',num2str(clumps),styleparam,' -l ',num2str(colRowName),' -L 0'];
    %disp(command)
    system(command);
    result = csvread('temp_out.csv', 1, 0);
    
end;

