#include "func.h"
#include "util.h"

#include "omp.h"

const int NTHREADS = 8;

void func0(double *weights, double *arrayX, double *arrayY, int xr, int yr, int n)
{
   // printf("Entered func0\n");
   int i;

   #pragma omp parallel for private(i) num_threads(NTHREADS)
   for(i = 0; i < n; i++){
      weights[i] = 1/((double)(n));
      arrayX[i] = xr;
      arrayY[i] = yr;
   }

   // printf("Exiting func0\n");

}

void func1(int *seed, int *array, double *arrayX, double *arrayY,
         double *probability, double *objxy, int *index,
         int Ones, int iter, int X, int Y, int Z, int n)
{
   // printf("Entered func1\n");

   int i, j;
      int index_X, index_Y;
   int max_size = X*Y*Z;

      for(i = 0; i < n; i++){
         arrayX[i] += 1 + 5*rand2(seed, i);
         arrayY[i] += -2 + 2*rand2(seed, i);
      }

      for(i = 0; i<n; i++){
         for(j = 0; j < Ones; j++){
            index_X = round(arrayX[i]) + objxy[j*2 + 1];
            index_Y = round(arrayY[i]) + objxy[j*2];
            index[i*Ones + j] = fabs(index_X*Y*Z + index_Y*Z + iter);
            if(index[i*Ones + j] >= max_size)
               index[i*Ones + j] = 0;
         }
         probability[i] = 0;

         for(j = 0; j < Ones; j++) {
            probability[i] += (pow((array[index[i*Ones + j]] - 100),2) -
                          pow((array[index[i*Ones + j]]-228),2))/50.0;
         }
         probability[i] = probability[i]/((double) Ones);
      }

   // printf("Exiting func1\n");

}

void func2(double *weights, double *probability, int n)
{
   // printf("Entered func2\n");   

   int i;
   double sumWeights=0;

   #pragma omp parallel for private(i) num_threads(NTHREADS)
   for(i = 0; i < n; i++)
         weights[i] = weights[i] * exp(probability[i]);

   #pragma omp parallel for private(i) num_threads(NTHREADS) reduction(+:sumWeights)
   for(i = 0; i < n; i++)
      sumWeights += weights[i];

   #pragma omp parallel for private(i) num_threads(NTHREADS)
   for(i = 0; i < n; i++)
         weights[i] = weights[i]/sumWeights;

   // printf("Exiting func2\n");

}

void func3(double *arrayX, double *arrayY, double *weights, double *x_e, double *y_e, int n)
{
   // printf("Entered func3\n");

   double estimate_x=0.0;
   double estimate_y=0.0;
    int i;

   for(i = 0; i < n; i++){
         estimate_x += arrayX[i] * weights[i];
         estimate_y += arrayY[i] * weights[i];
      }

   *x_e = estimate_x;
   *y_e = estimate_y;

   // printf("Exiting func3\n");


}

void func4(double *u, double u1, int n)
{
   // printf("Entered func4\n");

   int i;

   for(i = 0; i < n; i++){
         u[i] = u1 + i/((double)(n));
      }

   // printf("Exiting func4\n");


}

void func5(double *x_j, double *y_j, double *arrayX, double *arrayY, double *weights, double *cfd, double *u, int n)
{
   // printf("Entered func5\n");

   int i, j;

   for(j = 0; j < n; j++){
         //i = findIndex(cfd, n, u[j]);
         i = findIndexBin(cfd, 0, n, u[j]);
         if(i == -1)
            i = n-1;
         x_j[j] = arrayX[i];
         y_j[j] = arrayY[i];

      }

   for(i = 0; i < n; i++){
      arrayX[i] = x_j[i];
      arrayY[i] = y_j[i];
      weights[i] = 1/((double)(n));
   }

   // printf("Exiting func5\n");

}