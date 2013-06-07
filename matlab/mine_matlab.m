% 帮助文件
% 示例1：
% result=mine_matlab('Spellman.csv',0.6,15,0,'m',1);

% 示例2：
% result=mine_matlab('Spellman.csv',0.6,15,0,'p',[3 4]);

% 示例3：
% result=mine_matlab('Spellman.csv',0.6,15,0,'b',10);

% 示例4：
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
                disp('Variable ',num2str(params(1)),' V.S ',num2str(params(2)))
            else
                return;
            end;
        otherwise
            disp('Unknown method.')
    end
    if isunix
        z='/';
    else
        z='\';
    end;
    binpath='';
    if ismac
        binpath=[pwd,z,'RapidMic-macos'];
    end;
    if ispc
        binpath=[pwd,z,'RapidMic.exe'];
    end;
        
    command =[ binpath,' -i "',inputcsvfile,'" -o "temp_out.csv"',' -a ',num2str(alpha),' -c ',num2str(clumps),styleparam,' -l ',num2str(colRowName),' -L 0'];
    disp(command)
    system(command);
    result = csvread('temp_out.csv', 1, 0);
    
end;

