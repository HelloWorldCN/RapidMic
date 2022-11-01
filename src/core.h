#include <stdlib.h>

void PrintArrayd(double *x, int n);
void PrintArrayi(int *x, int n);
void PrintMatrixd(double **x, int m, int n);
void PrintMatrixi(int **x, int m, int n);
/* sort
 * Sort x and idx (of length n) according to x.
 */
void sort(double *x, int *idx, int n);

/* EquipartitionYAxis
 * Returns the map Q: D -> {0, ...,y-1}. See Algorithm 3 in SOM.
 *
 * Parameters
 *   Dy (IN): y-data sorted in increasing order
 *   n (IN): length of Dy
 *   y (IN): an integer greater than 1
 *   Qm (OUT): the map Q. Qm must be a preallocated vector
 *       of size n.
 * Return
 *   q : the real value of y. It can be < y.
 */
int EquipartitionYAxis(double *Dy, int n, int y, int *Qm);

/* GetSuperclumpsPartition
 * Returns the map P: D -> {0, ...,k-1}.
 *
 * Parameters
 *   Dx (IN) : x-data sorted in increasing order
 *   n (IN) : length of Dx
 *   Qm (IN) : the map Q computed by EquipartitionYAxis sorted
 *       in increasing order by Dx-values.
 *   k_hat (IN) : maximum number of clumps
 *       Pm (IN): the map P. Pm must be a preallocated vector
 *       of size n.
 * Return
 *   k : number of clumps in Pm.
 */
int GetSuperclumpsPartition(double *Dx, int n, int *Qm, int k_hat);

/* ApproxOptimizeXAxis
 * Returns the map P: D -> {0, ...,k-1}. See Algorithm 2 in SOM.
 *
 * Parameters
 *   Dx (IN) : x-data sorted in increasing order by Dx-values
 *   Dy (IN) : y-data sorted in increasing order by Dx-values
 *   n (IN) : length of Dx and Dy
 *   Qm (IN) : the map Q computed by EquipartitionYAxis sorted
 *       in increasing order by Dx-values.
 *   q (IN) : number of clumps in Qm
 *   Pm (IN) : the map P computed by GetSuperclumpsPartition
 *       sorted in increasing order by Dx-values.
 *   p (IN) : number of clumps in Pm
 *   x (IN) : grid size on x-values
 *   I (OUT) : the normalized mutual information vector. It
 *       will contain I_{k,2}, ..., I_{k, x}. I must be a
 *       preallocated array of dimension x-1.
 */
void ApproxOptimizeXAxis(double *Dx, double *Dy, int n, int *Qm, int q, int *Pm,
                         int p, int x, double *I);
