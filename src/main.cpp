//
//  main.cpp
//  MineMac
//
//  Created by tang on 13-1-8.
//  Copyright (c) 2013年 tang. All rights reserved.
//

#ifdef WIN32
#pragma comment(lib, "..\\lthread-win\\lib\\x86\\pthreadVC2.lib")
#include "..\libmine\cppmine.h"
#else
#include "../libmine/cppmine.h"
#endif
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <iostream>


#include <pthread.h>
#include "stringenc.h"





int main (int argc, char **argv)
{
	


	double PI;
	int i, n,m;
	//double *x, *y;

	int j;
	

	PI = 3.14159265;
	MINE *mine=new MINE();
	mine->run(argc,argv);
    mine->exportResult();


	/* build the problem */
	n = 17;
	m=800;
	double x[17] = {0,0,0,1,2,34,23,56,89,39,34,87,12,90,23,12,56};
	double y[17] = {78,0,0,23,56,1,65,96,6,2,23,3,13,5,57,2,8};
	double **inData=new double *[m];
	for (i=0;i<m;i++)
	{
		inData[i]=x;
	}
	//double *x = new double[n];
	//double *y = new double[n];
	//     for (i=0; i<n; i++)
	//     {
	//         /* build x = [0, 0.001, ..., 1] */
	//         x[i] = (double) i / (double) (n-1);
	// 		//x[i]=i;
	//
	//         /* build y = sin(10 * pi * x) + x */
	//         y[i] = sin(10 * PI * x[i]) + x[i];
	// 		//y[i]=-x[i];
	//
	//     }

	/* compute score */
	//mine->OnePairsAnalysis(x, y, n);

	/* print mine statistics */
	/*cout << "MIC: " << mine->get_mic() << "\n";
	cout << "MAS: " << mine->get_mas() << "\n";
	cout << "MEV: " << mine->get_mev() << "\n";
	cout << "MCN: " << mine->get_mcn() << "\n";*/

	//mine->AllPairsAnalysis((double **)inData, m,n);
	//mine->AllPairsAnalysis((double **)inData, m,n);
	//mine->TwoSetsAnalysis(inData, m, n, inData, m, n);
	/* delete the mine object */
	delete mine;
	delete []inData;
	//test();

	/* free the problem */
	//delete [] x;
	//delete [] y;
	//char c=getchar();
	//cout<< "tadfasdfds";

	return 0;
}