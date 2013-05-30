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


#define NUM_THREADS     5


using namespace std;
pthread_mutex_t count_lock= PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t count_nonzero=PTHREAD_COND_INITIALIZER;
int count=0;



void decrement_count()
{
	/*
	pthread_mutex_lock(&count_lock);
	while (count == 0)
	pthread_cond_wait(&count_nonzero, &count_lock);
	count = count - 1;
	pthread_mutex_unlock(&count_lock);
	*/
}


void* increment_count(void *threadid)
{
	/*long tid;
	tid = (long)threadid;
	int i;
	for (i=0;i<10000;i++)
	{
	printf("Hello World! It's me, thread #%ld!\n", tid);
	}
	pthread_mutex_lock(&count_lock);
	count = count + 1;
	if (count==5)
	{
	pthread_cond_signal(&count_nonzero);
	}
	pthread_mutex_unlock(&count_lock);
	pthread_exit(NULL);*/
	return NULL;
}

void testcondition ()
{
	pthread_t threads[NUM_THREADS];
	int res;

	void *thread_result[NUM_THREADS];
	int rc;
	long t;
	for(t=0; t<NUM_THREADS; t++){
		//printf("In main: creating thread %ld\n", t);
		rc = pthread_create(&threads[t], NULL, increment_count, (void *)t);

		if (rc){
			//printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

	}
	pthread_mutex_lock(&count_lock);

	pthread_cond_wait(&count_nonzero, &count_lock);

	pthread_mutex_unlock(&count_lock);
	printf("end!\n");
	/* Last thing that main() should do */
	pthread_mutex_destroy(&count_lock);
	pthread_cond_destroy(&count_nonzero);
	pthread_exit(NULL);
}
void *PrintHello(void *threadid)
{
	long tid;
	tid = (long)threadid;
	int i=0;
	for (i=0;i<10;i++)
	{
		printf("Hello World! It's me, thread #%ld!\n", tid);
	}

	pthread_exit(NULL);
	return NULL;
}

void test ()
{
	pthread_t threads[NUM_THREADS];
	int res;

	void *thread_result[NUM_THREADS];
	int rc;
	long t;
	for(t=0; t<NUM_THREADS; t++){
		//printf("In main: creating thread %ld\n", t);
		rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);

		if (rc){
			//printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

	}
	for(t=0; t<NUM_THREADS; t++){
		//printf("In main: creating thread %ld\n", t);

		res = pthread_join(threads[t], &thread_result[t]); //pthread_join 阻塞执行的线程直到某线程结束
		if (res != 0)
		{
			perror("Thread join failed");
			exit(EXIT_FAILURE);
		}
	}
	/* Last thing that main() should do */
	pthread_exit(NULL);
}


int main (int argc, char **argv)
{
	


	double PI;
	int i, n,m;
	//double *x, *y;

	int j;
	string str2("PLAYER,Record_ID#,SALARY,ROOKIE,POS,G,PA,AB,H,1B,2B,3B,HR,TB,BB,UBB,IBB,HBP,SF,SH,ROE,RBI,LEADOFF_PA,DP,TP,WP,PB,END_GAME,SO,BK,INTERFERENCE,FC,TOB,OUT,SIT_DP,GIDP,PITCHES,BALLS,STRIKES,FB,GB,LD,POP,Batted Balls,PA_P,PA_C,PA_1B,PA_2B,PA_3B,PA_SS,PA_LF,PA_CF,PA_RF,PA_DH,PA_PH,PA_PR,G_P,G_C,G_1B,G_2B,G_3B,G_SS,G_LF,G_CF,G_RF,G_DH,G_PH,G_PR,AVG,OBP,SLG,D_ISO,TBP,BBr,UBBr,IBBR,SO/BB,ABR,HITR,B1R,B2R,B3R,HRr,HBPR,SFR,SHR,ROEr,SOr,OUTR,NSOR,RBIR,LEADOFFR,END_GAMER,DP%,FB%,GB%,LD%,POP%,SB,CS,PICKOFF,R,SB%,RUNR,D_BPF,PA%,D_MLVr,D_PMLVr,D_RPMLVr,D_VORPr,D_MLV,D_PMLV,D_RPMLV,D_VORP,D_NETDP,D_EqA,D_EqR,D_RAR,D_RAP,D_RARP,D_OUTS_EQ,PA_ROB,R1,R2,R3,R1_BI,R2_BI,R3_BI,ROB,OBI,R1BI%,R2BI%,R3BI%,OBI%");
	vector<string> str1=split(str2,",");
	//划分任务
	for (i=0;i<5;i++)
	{
		for (j=i+1;j<5;j++)
		{

			//cout<<(i*(5-1)-i*(i-1)/2+j-i-1)<<" ";
			//outArray[]
		}
	}

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