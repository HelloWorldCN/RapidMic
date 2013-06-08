rem cl.bat
cl.exe main.cpp core.c mine.c cppmine.cpp arg_parser.cc stringenc.cpp /c /I".\lthread-win\include"
link.exe /out:main_cl.exe *.obj /LIBPATH:".\lthread-win\lib\x86" pthreadVC2.lib
pause