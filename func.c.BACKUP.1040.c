#include "func.h"
#include "util.h"

<<<<<<< Updated upstream
=======
#include "omp.h"
>>>>>>> Stashed changes

void func0(double *weights, double *arrayX, double *arrayY, int xr, int yr, int n)
{
	int i;
   #pragma omp parallel
   {
      #pragma omp for
   	for(i = 0; i < n; i++){
   		weights[i] = 1/((double)(n));
   		arrayX[i] = xr;
   		arrayY[i] = yr;
   	}
   }
}

#define NEST_FUNC1

void func1(int *seed, int *array, double *arrayX, double *arrayY,
			double *probability, double *objxy, int *index,
			int Ones, int iter, int X, int Y, int Z, int n)
{
	int i, j;
   	int index_X, index_Y;
	int max_size = X*Y*Z;

   #pragma omp parallel
   {

      #pragma omp for
   	for(i = 0; i < n; i++){
   		arrayX[i] += 1 + 5*rand2(seed, i);
   		arrayY[i] += -2 + 2*rand2(seed, i);
   	}

      #ifdef NEST_FUNC1
      #pragma omp for collapse(2) //2-nested for loop
      #else
      #pragma omp for
      #endif
      // each thread handles an (i,j) tuple
      // #pragma omp for //create team



      	for(i = 0; i<n; i++){
            #ifndef NEST_FUNC1
            // pull out repetitive multiplication and indexing
            int iOnes = i*Ones;
            double arrxi = arrayX[i];
            double arryi = arrayY[i];
            #endif
            // #pragma omp for //create team for each team
      		for(j = 0; j < Ones; j++){
               // pull out repetitive multiplication and indexing
               // TODO: is keeping this in nested loop to allow for
               // #pragma omp for collapse(2) (but having
               // extra multiplication/indexing)
               // faster or slower?
               #ifdef NEST_FUNC1
               int iOnes = i*Ones;
               double arrxi = arrayX[i];
               double arryi = arrayY[i];
               #endif

               // pull out repetitive multiplication
               int ioj = iOnes + j;
               int j2 = j*2;
      			index_X = round(arrxi) + objxy[j2 + 1];
      			index_Y = round(arryi) + objxy[j2];
      			index[ioj] = fabs(index_X*Y*Z + index_Y*Z + iter);
      			if(index[ioj] >= max_size)
      				index[ioj] = 0;
               // remove excess indexing
               int iioj = index[ioj];
               // combine for loops and remove need to initialize probability[i]
               probability[i] = (pow((array[iioj] - 100),2) -
                             pow((array[iioj]-228),2))/50.0;
               // probability[i] = (pow((array[index[i*Ones + j]] - 100),2) -
               //               pow((array[index[i*Ones + j]]-228),2))/50.0;
      		    #ifdef NEST_FUNC1
                probability[i] = probability[i]/((double) Ones);
                #endif
            }
      		// probability[i] = 0;

      		// for(j = 0; j < Ones; j++) {
      		// 	probability[i] += (pow((array[index[i*Ones + j]] - 100),2) -
      		// 					  pow((array[index[i*Ones + j]]-228),2))/50.0;
      		// }
            #ifndef NEST_FUNC1
      		probability[i] = probability[i]/((double) Ones);
      	   #endif
         }

   
   }
}

void func2(double *weights, double *probability, int n)
{
	int i;
	double sumWeights=0;

   #pragma omp parallel for //combine for loops
	  for(i = 0; i < n; i++)
      {
   		weights[i] = weights[i] * exp(probability[i]);
         sumWeights += weights[i];
         weights[i] = weights[i]/sumWeights;
      }

 //   	for(i = 0; i < n; i++)
 //   		sumWeights += weights[i];

	// for(i = 0; i < n; i++)
 //   		weights[i] = weights[i]/sumWeights;
   
}

void func3(double *arrayX, double *arrayY, double *weights, double *x_e, double *y_e, int n)
{
	double estimate_x=0.0;
	double estimate_y=0.0;
    int i;

    #pragma omp parallel for
   	for(i = 0; i < n; i++){
      		estimate_x += arrayX[i] * weights[i];
      		estimate_y += arrayY[i] * weights[i];
      	}
    
	*x_e = estimate_x;
	*y_e = estimate_y;

}

void func4(double *u, double u1, int n)
{
	int i;
<<<<<<< Updated upstream
   double dn = (double)(n);

=======
>>>>>>> Stashed changes
   #pragma omp parallel for
	for(i = 0; i < n; i++){
   		u[i] = u1 + i/dn;
   	}
}

void func5(double *x_j, double *y_j, double *arrayX, double *arrayY, double *weights, double *cfd, double *u, int n)
{
	int i, j;

   int weight = 1/((double)(n));

   #pragma omp parallel
   {
      #pragma omp for
   	for(j = 0; j < n; j++){
      		//i = findIndex(cfd, n, u[j]);
      		i = findIndexBin(cfd, 0, n, u[j]);
      		if(i == -1)
      			i = n-1;
      		x_j[j] = arrayX[i];
      		y_j[j] = arrayY[i];

      	}
      
      // hoist out division

      #pragma omp simd // all simple loads
   	for(i = 0; i < n; i++){
   		arrayX[i] = x_j[i];
   		arrayY[i] = y_j[i];
   		// weights[i] = 1/((double)(n));
         weights[i] = weight;
      	}
   }
}
