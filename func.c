#include "func.h"
#include "util.h"

#include "omp.h" //to actually have things work!

const int NUM_THREADS = 4;

void func0(double *weights, double *arrayX, double *arrayY, int xr, int yr, int n)
{
   printf("Entering func0\n");
	int i;
   // hoist out division
   double weight = 1/((double)(n));
   // #pragma omp parallel for default(shared) private(i) num_threads(NUM_THREADS)
   #pragma omp parallel
   {
      #pragma omp for private(i)
   	for(i = 0; i < n; i++){
   		weights[i] = weight;
   		arrayX[i] = xr;
   		arrayY[i] = yr;
   	}
   }
   printf("Exiting func0\n");

}

void func1(int *seed, int *array, double *arrayX, double *arrayY,
			double *probability, double *objxy, int *index,
			int Ones, int iter, int X, int Y, int Z, int n)
{
   printf("Entering func1\n");

	int i, j;
   	int index_X, index_Y;
	int max_size = X*Y*Z;

      // #pragma omp parallel for
   // #pragma omp parallel
   // {
      #pragma omp for private(i)
   	for(i = 0; i < n; i++){
   		arrayX[i] += 1 + 5*rand2(seed, i);
   		arrayY[i] += -2 + 2*rand2(seed, i);
   	}
   // }
      // #pragma omp parallel for
   // #pragma omp parallel 
   // {
   	for(i = 0; i<n; i++){
         // #pragma omp single //results can be shared during nested loop
         // {
         double arrxi = arrayX[i];
         double arryi = arrayY[i];
         int iones = i*Ones;
         // }
         // #pragma omp for private(i,j)
         #pragma omp for private(j)
   		for(j = 0; j < Ones; j++){
            // do less math by saving answer early
            int ioj = iones + j;
   			index_X = round(arrayX[i]) + objxy[j*2 + 1];
   			index_Y = round(arrayY[i]) + objxy[j*2];
            index[ioj] = fabs(index_X*Y*Z + index_Y*Z + iter);
            if(index[ioj] >= max_size)
               index[ioj] = 0;
            probability[i] = (pow((array[index[ioj]] - 100),2) -
                          pow((array[index[ioj]]-228),2))/50.0;
   			// index[i*Ones + j] = fabs(index_X*Y*Z + index_Y*Z + iter);
   			// if(index[i*Ones + j] >= max_size)
   			// 	index[i*Ones + j] = 0;
      //       probability[i] = (pow((array[index[i*Ones + j]] - 100),2) -
      //                     pow((array[index[i*Ones + j]]-228),2))/50.0;
   		}
   		// probability[i] = 0;

   		// for(j = 0; j < Ones; j++) {
   		// 	probability[i] += (pow((array[index[i*Ones + j]] - 100),2) -
   		// 					  pow((array[index[i*Ones + j]]-228),2))/50.0;
   		// }
   		probability[i] = probability[i]/((double) Ones);
   	// }
   }
   printf("Exiting func1\n");

}

void func2(double *weights, double *probability, int n)
{
   printf("Entering func2\n");

	int i;
	double sumWeights=0;

   // combine loops, parallelize

   // #pragma omp parallel for
   #pragma omp parallel
   { 
      int acc = 0; //private for each thread
      #pragma omp for private(i)
      for (i = 0; i < n; i++)
      {
         weights[i] *= exp(probability[i]);
         acc += weights[i];
         // sumWeights += weights[i];
         // weights[i] /= sumWeights;
      }
      #pragma omp atomic
      sumWeights += acc;
   }

   #pragma omp simd private(i)
   for (i = 0; i < n; i++)
   {
      weights[i] /= sumWeights;
   }

	// for(i = 0; i < n; i++)
 //   		weights[i] = weights[i] * exp(probability[i]);

 //   	for(i = 0; i < n; i++)
 //   		sumWeights += weights[i];

	// for(i = 0; i < n; i++)
 //   		weights[i] = weights[i]/sumWeights;
   printf("Exiting func2\n");

}

void func3(double *arrayX, double *arrayY, double *weights, double *x_e, double *y_e, int n)
{
   printf("Entering func3\n");

	double estimate_x=0.0;
	double estimate_y=0.0;
   int i;

   #pragma omp parallel
   {
    // #pragma omp parallel for reduction(+:estimate_x, estimate_y) private(i) num_threads(NUM_THREADS)
      int accx = 0;
      int accy = 0;
      #pragma omp for private(i)
   	for(i = 0; i < n; i++){
      		accx += arrayX[i] * weights[i];
      		accy += arrayY[i] * weights[i];
      	}
      #pragma omp atomic
      estimate_x += accx;
      #pragma omp atomic
      estimate_y += accy;
   }

	*x_e = estimate_x;
	*y_e = estimate_y;

   printf("Exiting func3\n");


}

void func4(double *u, double u1, int n)
{
   printf("Entering func4\n");

	int i;

   #pragma omp parallel for private(i)
	for(i = 0; i < n; i++){
   		u[i] = u1 + i/((double)(n));
   	}
   printf("Exiting func4\n");

}

void func5(double *x_j, double *y_j, double *arrayX, double *arrayY, double *weights, double *cfd, double *u, int n)
{
   printf("Entering func5\n");

	// int i, j;

   // hoist out division
   // int i,j;
   int i,j;

   // // #pragma omp parallel for private(j) num_threads(NUM_THREADS) //should not be allowed to share i, so declare inside loop
   // #pragma omp parallel
   // {
   //    #pragma omp for private(i,j)
      // #pragma omp parallel for private(i,j)
      // #pragma omp parallel for private(j)
      for(j = 0; j < n; j++){
   	// for(j = 0; j < n; j++){
      		//i = findIndex(cfd, n, u[j]);
            i = findIndexBin(cfd, 0, n, u[j]);
            // int i = findIndexBin(cfd, 0, n, u[j]);

      		// i = findIndexBin(cfd, 0, n, u[j]);
      		if(i == -1)
      			i = n-1;
      		x_j[j] = arrayX[i];
      		y_j[j] = arrayY[i];

      	}
   printf("Halfway through func5\n");

      // int i;
   double weight = 1/((double)(n));

      // #pragma omp parallel for private(i) num_threads(NUM_THREADS)
   	// for(i = 0; i < n; i++){
      // #pragma omp for private(i)
      #pragma omp parallel for private(i)
      for(i = 0; i < n; i++){
   		arrayX[i] = x_j[i];
   		arrayY[i] = y_j[i];
   		weights[i] = weight;
   	}
   // }
   printf("Exiting func5\n");

}
