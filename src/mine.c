
/*
 Dr. Tang rewrite the core MIC implementation based on Davide Albanese
 The following code support rapidly computing MIC.
 We presented a new rapidly computing maximal information-based nonparametric
 exploration tool for statistical analysis of large-scale dataset.  By parallel
 processing of MIC algorithm, the algorithm can well analyze large-scale dataset
 and greatly reduce the coputing time. 下面的代码由西南交通大学 Dr唐 改写
 本代码采用多线程的方式快速计算MIC 值

 Southwest Jiaotong University

 */
/*
 This code is written by Davide Albanese <davide.albanese@gmail.com>.
 (C) 2012 Davide Albanese, (C) 2012 Fondazione Bruno Kessler.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mine.h"
#include "core.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(WIN32) || defined(_WIN32)
#pragma comment(lib, ".\\lthread-win\\lib\\x86\\pthreadVC2.lib")
#include ".\\lthread-win\\include\\pthread.h"
#else
#include <pthread.h>
#endif
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX_THREADS_NUMS 50
#define MINE_SUCCESS 1
#define MINE_FAIL 0

pthread_mutex_t thread_complete_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t res_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t all_thread_completed = PTHREAD_COND_INITIALIZER;
unsigned thread_complete_count = 0;

pthread_mutex_t batch_thread_complete_lock = PTHREAD_MUTEX_INITIALIZER;
// pthread_mutex_t batch_res_lock= PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t batch_all_thread_completed = PTHREAD_COND_INITIALIZER;
unsigned batch_thread_complete_count = 0;

int getOptimizedThreadNum(int dataNum, int flag) {
  // for single partion
  if (flag) {
    // linear fitting
    return MIN(MAX_THREADS_NUMS,
               (int)floor(0.002152803091055 * dataNum + 33.243420573011576));
  } else {
  }
  return 0;
}

void mine_problem_init(mine_problem *prob, mine_parameter *param) {
  int i, j;
  double B;
  int gx, gy, Gy_max;
  prob->score = (mine_score *)malloc(sizeof(mine_score));
  prob->score->m = 0;
  /* x-1 for each y */
  prob->score->p = NULL;
  prob->score->I = NULL;
  prob->Gy = NULL;

  B = MAX(pow(prob->n, param->alpha), 4);
  Gy_max = MAX((int)floor(B / 2.0), 2);

  prob->score->m = Gy_max - 1; //
  prob->score->p = (int *)malloc(prob->score->m * sizeof(int));
  prob->Gy = (int *)malloc(prob->score->m * sizeof(int));
  for (gy = 2; gy < Gy_max + 1; gy++) {
    gx = (int)floor(B / gy);
    if ((gx * gy) <= B) {
      prob->score->p[gy - 2] =
          gx -
          1; //最小划分为2格，如果gx为1则不用划分了，所以也就不用存储数据，如果为2则只有一个数据
      prob->Gy[gy - 2] = gy;
    }
  }

  prob->score->I = (double **)malloc(prob->score->m * sizeof(double *));

  /* sort and argsort of x and y */
  prob->x_x = (double *)malloc(prob->n * sizeof(double));
  prob->y_y = (double *)malloc(prob->n * sizeof(double));
  prob->idx_x = (int *)malloc(prob->n * sizeof(int));
  prob->idx_y = (int *)malloc(prob->n * sizeof(int));
  for (i = 0; i < prob->n; i++) {
    prob->x_x[i] = prob->x[i];
    prob->y_y[i] = prob->y[i];
    prob->idx_x[i] = i;
    prob->idx_y[i] = i;
  }
  //对原始序列值进行排序，同时记录排序后的原数据所处位置序列
  sort(prob->x_x, prob->idx_x, prob->n);
  sort(prob->y_y, prob->idx_y, prob->n);

  /* sort x by y-value and sort y by x-value */
  prob->x_y = (double *)malloc(
      prob->n * sizeof(double)); //对(xi,yi)的成对数据按照y的值进行排序后的x序列
  prob->y_x = (double *)malloc(
      prob->n * sizeof(double)); //对(xi,yi)的成对数据按照x的值进行排序后的y序列
  for (j = 0; j < prob->n; j++) {
    prob->x_y[j] = prob->x[prob->idx_y[j]];
    prob->y_x[j] = prob->y[prob->idx_x[j]];
  }
}

//************************************
// Method:    FixOnePartition
// 当x或者是y划分成固定的格数时，求另外一个数据划分成不同的格数时的信息
// FullName:  FixOnePartition
// Access:    public
// Returns:   void
// Qualifier:
// Parameter: mine_problem * prob
// Parameter: mine_parameter * param
// Parameter: int firstpartition
// Parameter: int i
// Parameter: double * It
//************************************
void FixOnePartition(mine_problem *prob, double c, int firstpartition,
                     int indexp, double *It, int *Qm, int *Qm_s) {
  int k;
  int j, q, p;

  /* x vs. y */
  //这里可以转换为多个线程同步进行，当y的划分个数固定的时候，通过动态规划方式确定不同的x划分的方式
  if (firstpartition == 0) {
    k = MAX((int)(c * (prob->score->p[indexp] + 1)), 1);
    /*等分为最多Gy[indexp]份，但是要保证相同值的数据在一个格子中，因此最终划分的个数可能小于Gy[indexp]*/
    q = EquipartitionYAxis(prob->y_y, prob->n, prob->Gy[indexp], Qm_s); //
    // 		PrintArrayi(idx_x,prob->n);
    // 		PrintArrayi(idx_y,prob->n);
    // 		PrintArrayi(Qm_s,prob->n);

    for (j = 0; j < prob->n; j++)
      Qm[prob->idx_y[j]] = Qm_s[j]; //获得y在原始排列情况下对应的划分表

    for (j = 0; j < prob->n; j++)
      Qm_s[j] = Qm
          [prob->idx_x
               [j]]; //获得x在原始数据对(x,y)约束条件下x排序后数据划分表,也可以理解成和x对位的y所在划分格子

    memcpy(Qm, Qm_s, sizeof(int) * prob->n);

    //在y划定的情况下重新确定x的分布，要保证x排序后的数据中值相同的数据要在一个格子中
    p = GetSuperclumpsPartition(prob->x_x, prob->n, Qm, k);
    /*PrintArrayi(score->p,score->m);*/
    ApproxOptimizeXAxis(prob->x_x, prob->y_x, prob->n, Qm_s, q, Qm, p,
                        prob->score->p[indexp] + 1, It);
  }

  if (firstpartition == 1) {
    k = MAX((int)(c * (prob->score->p[indexp] + 1)), 1);
    q = EquipartitionYAxis(prob->x_x, prob->n, prob->Gy[indexp], Qm_s);
    for (j = 0; j < prob->n; j++)
      Qm[prob->idx_x[j]] = Qm_s[j];

    for (j = 0; j < prob->n; j++)
      Qm_s[j] = Qm[prob->idx_y[j]];
    memcpy(Qm, Qm_s, sizeof(int) * prob->n);
    p = GetSuperclumpsPartition(prob->y_y, prob->n, Qm, k);

    ApproxOptimizeXAxis(prob->y_y, prob->x_y, prob->n, Qm_s, q, Qm, p,
                        prob->score->p[indexp] + 1, It);
  }
}

/*This function is the work thread definition and implementation for finding the
 highest mutual information attainable with multi different partion grids on
 variable sequence x and y*/
// pparam contains multi different partion grids
void *PartitionThread(void *pparam) {
  int i = 0, j = 0;
  threadparams *pthdparainf = (threadparams *)pparam;
  int *Qm, *Qm_s;
  double **subI;
  Qm = (int *)malloc(pthdparainf->prob->n * sizeof(int));
  Qm_s = (int *)malloc(pthdparainf->prob->n * sizeof(int));

  subI = (double **)malloc((pthdparainf->indexEnd - pthdparainf->indexBegin) *
                           sizeof(double *));
  for (i = pthdparainf->indexBegin; i < pthdparainf->indexEnd; i++) {
    subI[i - pthdparainf->indexBegin] =
        (double *)malloc((pthdparainf->prob->score->p[i]) * sizeof(double));
    FixOnePartition(pthdparainf->prob, pthdparainf->c,
                    pthdparainf->firstPartition, i,
                    subI[i - pthdparainf->indexBegin], Qm, Qm_s);
  }

  // avoid frequently mutually operation
  pthread_mutex_lock(&res_lock);
  for (i = pthdparainf->indexBegin; i < pthdparainf->indexEnd; i++) {
    if (pthdparainf->firstPartition == 0) {
      // compare i row
      for (j = 0; j < pthdparainf->prob->score->p[i]; j++)
        pthdparainf->prob->score->I[i][j] =
            MAX(subI[i - pthdparainf->indexBegin][j],
                pthdparainf->prob->score->I[i][j]);
    } else if (pthdparainf->firstPartition == 1) {
      // compare i column
      for (j = 0; j < pthdparainf->prob->score->p[i]; j++)
        pthdparainf->prob->score->I[j][i] =
            MAX(subI[i - pthdparainf->indexBegin][j],
                pthdparainf->prob->score->I[j][i]);
    }
    free(subI[i - pthdparainf->indexBegin]);
  }
  pthread_mutex_unlock(&res_lock);
  free(Qm);
  free(Qm_s);
  free(subI);

  pthread_mutex_lock(&thread_complete_lock);
  thread_complete_count--;
  if (thread_complete_count == 0) {
    pthread_cond_signal(&all_thread_completed);
  }
  pthread_mutex_unlock(&thread_complete_lock);
  pthread_exit(NULL);
  return NULL;
}

// compute one pair variables using multi-threads
mine_score *mine_compute_score_mt(mine_problem *prob, mine_parameter *param) {
  int i = 0, j = 0;
  threadparams *pthreadinfos;
  pthread_t *threads;
  int pageSize;
  int t, rc, numThreads;

  /* x vs. y */
  for (i = 0; i < prob->score->m; i++) {
    // pre-allocate space
    prob->score->I[i] = (double *)malloc((prob->score->p[i]) * sizeof(double));
    for (j = 0; j < prob->score->p[i]; j++) {
      prob->score->I[i][j] = 0.0f;
    }
  }
  numThreads = getOptimizedThreadNum(prob->n, 1);
  pthreadinfos = (threadparams *)malloc(sizeof(threadparams) * numThreads);
  threads = (pthread_t *)malloc(sizeof(pthread_t) * numThreads);
  pageSize = (2 * prob->score->m) / numThreads;
  thread_complete_count = numThreads;
  for (t = 0; t < numThreads; t++) {
    pthreadinfos[t].c = param->c;
    if (t < (numThreads / 2 - 1)) {
      pthreadinfos[t].firstPartition = 0;
      pthreadinfos[t].indexBegin = (t)*pageSize;
      pthreadinfos[t].indexEnd = (t + 1) * pageSize;
      pthreadinfos[t].prob = prob;
    }
    if (t == (numThreads / 2 - 1)) {
      pthreadinfos[t].firstPartition = 0;
      pthreadinfos[t].indexBegin = (t)*pageSize;
      pthreadinfos[t].indexEnd = prob->score->m;
      pthreadinfos[t].prob = prob;
    }
    if (t > (numThreads / 2 - 1) && t < (numThreads - 1)) {
      pthreadinfos[t].firstPartition = 1;
      pthreadinfos[t].indexBegin = (t - numThreads / 2) * pageSize;
      pthreadinfos[t].indexEnd = (t + 1 - numThreads / 2) * pageSize;
      pthreadinfos[t].prob = prob;
    }
    if (t == (numThreads - 1)) {
      pthreadinfos[t].firstPartition = 1;
      pthreadinfos[t].indexBegin = (t - numThreads / 2) * pageSize;
      pthreadinfos[t].indexEnd = prob->score->m;
      pthreadinfos[t].prob = prob;
    }
    rc = pthread_create(&threads[t], NULL, PartitionThread,
                        (void *)&pthreadinfos[t]);

    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  pthread_mutex_lock(&thread_complete_lock);
  pthread_cond_wait(&all_thread_completed, &thread_complete_lock);
  pthread_mutex_unlock(&thread_complete_lock);

  // pthread_exit(NULL);
  return prob->score;
}

//
// This function is the work thread definition and implementation for batch
// computing multi pairs variables pparam contains multi pais variables
// informations
void *batchComputeScoreThread(void *pparam) {
  int i = 0, j = 0, k = 0;
  mine_problem prob;
  mine_score score;
  mine_parameter *param;
  batchThreadparams *pthdparainf = (batchThreadparams *)pparam;
  int *Qm, *Qm_s;
  double *It;
  param = pthdparainf->param;

  // prob=(mine_problem *)malloc(sizeof(mine_problem));
  prob.n = pthdparainf->n;
  /* x-1 for each y */

  prob.score = &score;
  score.m = pthdparainf->scorem;
  score.p = pthdparainf->scorep;

  prob.x_x = (double *)malloc(prob.n * sizeof(double));
  prob.y_y = (double *)malloc(prob.n * sizeof(double));
  prob.idx_x = (int *)malloc(prob.n * sizeof(int));
  prob.idx_y = (int *)malloc(prob.n * sizeof(int));
  prob.x_y = (double *)malloc(
      prob.n * sizeof(double)); //对(xi,yi)的成对数据按照y的值进行排序后的x序列
  prob.y_x = (double *)malloc(
      prob.n * sizeof(double)); //对(xi,yi)的成对数据按照x的值进行排序后的y序列
  Qm = (int *)malloc(prob.n * sizeof(int)); //原始数据排列顺序划分表
  Qm_s =
      (int *)malloc(prob.n * sizeof(int)); //按照数据排序后顺序保存的数据划分表
  score.I = (double **)malloc(score.m * sizeof(double *));
  prob.Gy = pthdparainf->Gy;
  for (i = 0; i < pthdparainf->parisLength; i++) {
    prob.x = pthdparainf->inData[pthdparainf->paris[i].var1];
    prob.y = pthdparainf->inData[pthdparainf->paris[i].var2];

    /* sort and argsort of x and y */

    for (j = 0; j < prob.n; j++) {
      if (i == 0 ||
          (i != 0 && pthdparainf->paris[i].var1 !=
                         pthdparainf->paris[i - 1]
                             .var1)) { //和上一轮不相同则需要重新赋值并排序
        prob.x_x[j] = prob.x[j];
        prob.idx_x[j] = j;
      }

      if (i == 0 ||
          (i != 0 && pthdparainf->paris[i].var2 !=
                         pthdparainf->paris[i - 1]
                             .var2)) { //和上一轮不相同则需要重新赋值并排序
        prob.y_y[j] = prob.y[j];
        prob.idx_y[j] = j;
      }
    }
    /* sort x by y-value and sort y by x-value */
    //对原始序列值进行排序，同时记录排序后的原数据所处位置序列
    if (i == 0 || (i != 0 && pthdparainf->paris[i].var1 !=
                                 pthdparainf->paris[i - 1].var1))
      sort(prob.x_x, prob.idx_x, prob.n);
    if (i == 0 || (i != 0 && pthdparainf->paris[i].var2 !=
                                 pthdparainf->paris[i - 1].var2))
      sort(prob.y_y, prob.idx_y, prob.n);

    for (j = 0; j < prob.n; j++) {
      prob.x_y[j] = prob.x[prob.idx_y[j]];
      prob.y_x[j] = prob.y[prob.idx_x[j]];
    }

    /* x vs. y */
    //这里可以转换为多个线程同步进行，当y的划分个数固定的时候，通过动态规划方式确定不同的x划分的方式
    for (j = 0; j < score.m; j++) {
      score.I[j] = (double *)malloc((score.p[j]) * sizeof(double));
      FixOnePartition(&prob, param->c, 0, j, score.I[j], Qm, Qm_s);
    }

    /* y vs. x */
    for (k = 0; k < score.m; k++) {
      It = (double *)malloc((score.p[k]) * sizeof(double));
      FixOnePartition(&prob, param->c, 1, k, It, Qm, Qm_s);
      for (j = 0; j < score.p[k]; j++)
        score.I[j][k] = MAX(It[j], score.I[j][k]);

      free(It);
    }

    switch (pthdparainf->styleType) {
    case AllParis:
      get_result_score(&score, &pthdparainf->outArray[i]);
      break;
    case TwoSets:
      get_result_score(&score, &pthdparainf->outArray[i]);
      break;
    case MasterVariable:
      get_result_score(&score,
                       &pthdparainf->outArray[pthdparainf->paris[i].var1]);
      break;
    default:
      break;
    }

    for (j = 0; j < score.m; j++) {
      free(score.I[j]);
    }
  }
  free(Qm);
  free(Qm_s);
  free(prob.x_x);
  free(prob.y_y);
  free(prob.idx_x);
  free(prob.idx_y);
  free(prob.x_y);
  free(prob.y_x);
  free(score.I);

  pthread_mutex_lock(&batch_thread_complete_lock);
  batch_thread_complete_count = batch_thread_complete_count - 1;
  if (batch_thread_complete_count == 0) {
    pthread_cond_signal(&batch_all_thread_completed);
  }
  pthread_mutex_unlock(&batch_thread_complete_lock);

  pthread_exit((void *)0);
  return NULL;
}

/* Computes the maximum normalized mutual information scores
 * and returns a mine_score structure.
 */
mine_score *mine_compute_score(mine_problem *prob, mine_parameter *param) {
  int i = 0, j = 0;
  int *Qm, *Qm_s;
  double *It;
  Qm = (int *)malloc(prob->n * sizeof(int)); //原始数据排列顺序划分表
  Qm_s =
      (int *)malloc(prob->n * sizeof(int)); //按照数据排序后顺序保存的数据划分表

  /* x vs. y */
  //这里可以转换为多个线程同步进行，当y的划分个数固定的时候，通过动态规划方式确定不同的x划分的方式
  for (i = 0; i < prob->score->m; i++) {
    prob->score->I[i] = (double *)malloc((prob->score->p[i]) * sizeof(double));
    FixOnePartition(prob, param->c, 0, i, prob->score->I[i], Qm, Qm_s);
  }

  /* y vs. x */
  for (i = 0; i < prob->score->m; i++) {
    It = (double *)malloc((prob->score->p[i]) * sizeof(double));
    FixOnePartition(prob, param->c, 1, i, It, Qm, Qm_s);
    for (j = 0; j < prob->score->p[i]; j++)
      prob->score->I[j][i] = MAX(It[j], prob->score->I[j][i]);
    free(It);
  }

  free(Qm);
  free(Qm_s);

  return prob->score;
}

/* This function checks the parameters. It should be called
 * before calling mine_compute_score(). It returns NULL if
 * the parameters are feasible, otherwise an error message is returned.
 */
char *check_parameter(mine_parameter *param) {
  if ((param->alpha <= 0.0) || (param->alpha > 1.0))
    return "alpha must be in (0, 1.0]";

  if (param->c <= 0.0)
    return "c must be > 0.0";

  return NULL;
}

/*
 This function get mic,mev,mas,mcn from mine_score
 */
void get_result_score(mine_score *score, mine_result_score *result) {
  int i, j, b;
  double score_max_mic, s;
  double score_max_mas, score_max_mev, score_max_mcn;
  int b_max = 4;
  score_max_mcn = -1.0;
  score_max_mic = score_max_mev = score_max_mas = 0.0;

  if (score != NULL) {
    for (i = 0; i < score->m; i++)
      for (j = 0; j < score->p[i]; j++) {
        if (score->I[i][j] > score_max_mic)
          score_max_mic = score->I[i][j];
        s = fabs(score->I[i][j] - score->I[j][i]);
        if (s > score_max_mas)
          score_max_mas = s;
        if (((j == 0) || (i == 0)) && score->I[i][j] > score_max_mev)
          score_max_mev = score->I[i][j];
        b = (i + 2) * (j + 2);
        if ((score->I[i][j] > score_max_mcn) ||
            ((score->I[i][j] == score_max_mcn) && (b < b_max))) {
          score_max_mcn = score->I[i][j];
          b_max = b;
        }
      }
  }
  result->mic = score_max_mic;
  result->mas = score_max_mas;
  result->mev = score_max_mev;
  result->mcn = log(b_max) / log(2.0);
}

/* This function frees the memory used by a mine_score and
 *  destroys the score structure.
 */
void mine_free_score(mine_score **score) {
  int i;
  mine_score *score_ptr = *score;

  if (score_ptr != NULL) {
    if (score_ptr->m != 0) {
      free(score_ptr->p);
      for (i = 0; i < score_ptr->m; i++)
        free(score_ptr->I[i]);
      free(score_ptr->I);
    }

    free(score_ptr);
    score_ptr = NULL;
  }
}

// This function create multi work threads for analyzing input data rely on task
// type.
int createBatchComputeThread(mine_parameter *param, double **inDataSet, int m,
                             int n, int masterOrBetweenVariableid,
                             mine_result_score *outArray, int outLen,
                             ANALYSISSTYLES styleType) {
  int i, j, t, k, rc;
  int numThreads;
  //划分任务
  batchThreadparams *pthreadinfos;
  pthread_t *threads;
  void **thread_result;
  int pageSize;
  varPairs *pVarPairs;
  int index;
  double B;
  int gx, gy, Gy_max, forjbegin, forjend, foriend;
  int *p, *Gy;

  B = MAX(pow(n, param->alpha), 4);
  Gy_max = MAX((int)floor(B / 2.0), 2);
  p = (int *)malloc((Gy_max - 1) * sizeof(int));
  Gy = (int *)malloc((Gy_max - 1) * sizeof(int));
  for (gy = 2; gy < Gy_max + 1; gy++) {
    gx = (int)floor(B / gy);
    if ((gx * gy) <= B) {
      p[gy - 2] =
          gx -
          1; //最小划分为2格，如果gx为1则不用划分了，所以也就不用存储数据，如果为2则只有一个数据
      Gy[gy - 2] = gy;
    }
  }
  if (outLen > 10000) {
    numThreads = 50;
  } else if (outLen > 5000) {
    numThreads = 25;
  } else if (outLen > 2500) {
    numThreads = 10;
  } else if (outLen > 1000) {
    numThreads = 5;
  } else if (outLen > 100) {
    numThreads = 2;
  } else
    numThreads = 1;

  batch_thread_complete_count = numThreads;
  pthreadinfos =
      (batchThreadparams *)malloc(sizeof(batchThreadparams) * numThreads);
  threads = (pthread_t *)malloc(sizeof(pthread_t) * numThreads);
  thread_result = (void **)malloc(sizeof(void *) * numThreads);

  if (styleType == MasterVariable)
    outLen--;

  pageSize = (int)floor((double)outLen / (double)numThreads);

  t = 0;
  k = 0;
  pVarPairs = (varPairs *)malloc(sizeof(varPairs) * pageSize);
  switch (styleType) {
  case AllParis:
    foriend = forjend = m;
    break;
  case TwoSets:
    forjbegin = masterOrBetweenVariableid;
    forjend = m;
    foriend = masterOrBetweenVariableid;
    break;
  case MasterVariable:
    foriend = m;
    forjend = 1;
    forjbegin = 0;
    break;
  default:
    break;
  }
  for (i = 0; i < foriend; i++) {
    if (styleType == AllParis) {
      forjbegin = i + 1;
    }
    if (styleType == MasterVariable && i == masterOrBetweenVariableid) {
      continue;
    }
    for (j = forjbegin; j < forjend; j++) {
      if (t < pageSize) {
        pVarPairs[t].var1 = i;
        if (styleType == MasterVariable) {
          pVarPairs[t].var2 = masterOrBetweenVariableid;
        } else
          pVarPairs[t].var2 = j;
        t++;
      }
      if (t == pageSize) {
        // pVarPairs=(varPairs *)malloc(sizeof(varPairs)*pageSize);
        switch (styleType) {
        case AllParis:
          pthreadinfos[k].styleType = AllParis;
          pthreadinfos[k].inData = inDataSet;
          index = i * (m - 1) - i * (i - 1) / 2 + j - i - 1 - pageSize + 1;
          pthreadinfos[k].outArray = &outArray[index];
          break;
        case TwoSets:
          pthreadinfos[k].styleType = TwoSets;
          pthreadinfos[k].inData = inDataSet;
          if (k < (numThreads - 1)) {
            pthreadinfos[k].outArray = &outArray[k * pageSize];
          } else {
            pthreadinfos[k].outArray = &outArray[outLen - pageSize];
          }

          break;
        case MasterVariable:
          pthreadinfos[k].styleType = MasterVariable;
          pthreadinfos[k].inData = inDataSet;
          pthreadinfos[k].outArray = outArray;
          break;
        default:
          break;
        }
        pthreadinfos[k].Gy = Gy;
        pthreadinfos[k].n = n;
        pthreadinfos[k].param = param;
        pthreadinfos[k].paris = pVarPairs;
        pthreadinfos[k].parisLength = t;
        pthreadinfos[k].scorem = Gy_max - 1;
        pthreadinfos[k].scorep = p;
        rc = pthread_create(&threads[k], NULL, batchComputeScoreThread,
                            (void *)&pthreadinfos[k]);
        // sleep(1);
        if (0 != rc)
          exit(-1);

        k++;
        t = 0;
        if (k < (numThreads - 1)) {
          pVarPairs = (varPairs *)malloc(sizeof(varPairs) * pageSize);
        } else if (k == (numThreads - 1)) {
          pageSize = outLen - k * pageSize;
          pVarPairs = (varPairs *)malloc(sizeof(varPairs) * pageSize);
        }
      }
    }
  }
  pthread_mutex_lock(&batch_thread_complete_lock);
  pthread_cond_wait(&batch_all_thread_completed, &batch_thread_complete_lock);
  pthread_mutex_unlock(&batch_thread_complete_lock);

  //回收资源

  free(Gy);
  free(p);
  free(threads);

  for (i = 0; i < numThreads; i++) {
    free(pthreadinfos[i].paris);
  }
  free(pthreadinfos);
  free(thread_result);
  // pthread_exit(NULL);
  return MINE_SUCCESS;
}

int mine_allPairs_analysis(mine_parameter *param, double **inData, int m, int n,
                           mine_result_score *outArray, int outLen) {
  return createBatchComputeThread(param, inData, m, n, 0, outArray, outLen,
                                  AllParis);
}

//
int mine_twoSetsAnalysis(mine_parameter *param, double **inDataSet, int m,
                         int n, int betweenid, mine_result_score *outArray,
                         int outLen) {
  return createBatchComputeThread(param, inDataSet, m, n, betweenid, outArray,
                                  outLen, TwoSets);
}

int mine_masterVariableAnalysis(mine_parameter *param, double **inData, int m,
                                int n, int masterid,
                                mine_result_score *outArray, int outLen) {
  return createBatchComputeThread(param, inData, m, n, masterid, outArray,
                                  outLen, MasterVariable);
}

int mine_onePair_analysis(mine_parameter *param, double *x, double *y, int n,
                          mine_result_score *outArray) {
  mine_problem *prob = (mine_problem *)malloc(sizeof(mine_problem));
  prob->x = x;
  prob->y = y;
  prob->n = n;
  if (prob->n >= 100) // Avoid expensive multi threads calls.
  {
    mine_problem_init(prob, param);
    mine_compute_score_mt(prob, param);
  } else {
    mine_problem_init(prob, param);
    mine_compute_score(prob, param);
  }

  get_result_score(prob->score, outArray);
  mine_free_score(&prob->score);
  free(prob->x_x);
  free(prob->y_y);
  free(prob->idx_x);
  free(prob->idx_y);
  free(prob->x_y);
  free(prob->y_x);
  free(prob);
  return MINE_SUCCESS;
}