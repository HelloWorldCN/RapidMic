//
//  main.cpp
//  MineMac
//
//  Created by tang on 13-1-8.
//  Copyright (c) 2013年 tang. All rights reserved.
//

#include "cppmine.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#if defined(WIN32) || defined(_WIN32)
#include <time.h>
#include <windows.h>
#else
#include <sys/time.h>
#endif

int main(int argc, char **argv) {

  MINE *mine = new MINE();
  mine->run(argc, argv);
  delete mine;

  return 0;
}