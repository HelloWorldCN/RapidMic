//
//  main.cpp
//  MineMac
//
//  Created by tang on 13-1-8.
//  Copyright (c) 2013年 tang. All rights reserved.
//

#include "cppmine.h"
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <iostream>

#if defined(WIN32) || defined(_WIN32)
#   include <windows.h>
#include <time.h>
#else
#   include <sys/time.h>
#endif



int main (int argc, char **argv)
{

	MINE *mine=new MINE();
	mine->run(argc,argv);
	delete mine;

	return 0;
}