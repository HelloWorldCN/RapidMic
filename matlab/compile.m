command='';
if isunix
    command='g++ ../src/main.cpp ../src/core.c ../src/mine.c ../src/cppmine.cpp ../src/arg_parser.cc ../src/stringenc.cpp -o RapidMic -pthread -O3 -D_MATLAB';
end;



%disp(command)
system(command);