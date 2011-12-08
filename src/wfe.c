#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"


// To do: change C.it into double

void Index(int *index, int *uniq_index, int *len_u_index, int *len_data, int *result){
  
  // initializing index_result
  memset(result, 0, sizeof(int)*(*len_data));

  int i;
  for (i = 0 ; i < *len_u_index; i++) {
    
    int k;
    for (k = 0 ; k < *len_data; k++) {
      if ( index[k] == uniq_index[i] ) {
	result[k] = i+1;
	//Rprintf("allocated index: %d\n", result[k]);
      }
    }

  }

}

void Vectorize(double *Wvec, int *nrow, int *ncol, int *time_index, int *dyad_index, int *n_obs, double *results) {

  int i, j, itemp;
  double **W = doubleMatrix(*nrow, *ncol);

  itemp = 0;
  for (j = 0; j < *ncol; j++)
    for (i = 0; i < *nrow; i++)
      W[i][j] = Wvec[itemp++];
  
  for (i = 0; i < *n_obs; i++) {
    results[i] = W[time_index[i]-1][dyad_index[i]-1];
  }
  
  FreeMatrix(W, *nrow);
}



void Transform(double *y, int *n, int *treat, double *pscore, double *ytrans) {

  int i, sumTreat;
  double psDenom0, psDenom1;

  sumTreat = 0; psDenom1 = 0; psDenom0 = 0;
  for (i = 0; i < *n; i++) {
    sumTreat += treat[i];
    if (treat[i] == 1) {
      psDenom1 += (1/pscore[i]);
    } else { 
      psDenom0 += (1/(1-pscore[i]));
    }
  }

  for (i = 0; i < *n; i++) {
    if (treat[i] == 1) {
      ytrans[i] = y[i]*sumTreat / (pscore[i] * psDenom1);
    } else {
      ytrans[i] = y[i]*(*n - sumTreat) / ((1 - pscore[i]) * psDenom0);      
    }
  }

}


// Generating time index when time index is missing
void GenTime(int* unit_index, 
    int* len_data, /* total number of rows in the data */
    int* len_u_index,
    double* time_index) {
  
  int i;
  int k;
  int count;
  for (i = 0 ; i < *len_u_index; i++) {
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( unit_index[k] == (i+1) ) {
        count ++;
	time_index[k] = count;
      }
    }
  } 
}




// Demeaning numerical vector by index
void Demean(double *var,
    int *index, 
    int *len_index, /* unique number of unit index, e.g. number of countries*/
    int *len_data, /* total number of rows of data */
    double *demean) {


  double count;
  double sum_var;
  double *mean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    sum_var = 0;
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        count ++;
        sum_var += var[k];
      }
    }
    mean[i] = sum_var / count;
    // Rprintf("1] Wmean_i: %f\n", Wmean[i]);
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        demean[k] = var[k] - mean[i];
        //Rprintf("2] dmean_k: %f\n", demean[k]);
      }
    }
  } 
}
    





// Weighted Demeaning numerical vector by index


void WDemean(double *var, double *weight,
    int *index, 
    int *len_index, /* unique number of unit index, e.g. number of countries*/
    int *len_data, /* total number of rows of data */
    double *Wdemean) {


  double sum_weight;
  double Wsum_var;
  double *Wmean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    Wsum_var = 0;
    sum_weight = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        sum_weight += weight[k];
        Wsum_var += weight[k] * var[k];
      }
    }
    Wmean[i] = Wsum_var / sum_weight;
    // Rprintf("1] Wmean_i: %f\n", Wmean[i]);
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        Wdemean[k] = var[k] - Wmean[i];
        // Rprintf("2] Wmean_k: %f\n", Wdemean[k]);
      }
    }
  } 
}
    


// Sqrt(W_it) * Weighted Demean_it


void WWDemean(double *var, double *weight,
    int *index, 
    int *len_index, /* unique number of unit index, e.g. number of countries*/
    int *len_data, /* total number of rows of data */
    double *WWdemean) {


  double sum_weight;
  double Wsum_var;
  double *Wmean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    Wsum_var = 0;
    sum_weight = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        sum_weight += weight[k];
        Wsum_var += weight[k] * var[k];
      }
    }
    Wmean[i] = Wsum_var / sum_weight;
    // Rprintf("1] Wmean_i: %f\n", Wmean[i]);
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
        WWdemean[k] = sqrt(weight[k]) * (var[k] - Wmean[i]);
        // Rprintf("3] WWmean_k: %f\n", WWdemean[k]);
      }
    }
  } 
}
    



int is_time_index_exist(int* u_i, int* t_i, int i, int j, int size) {
  int iter = 0;
  int exist = 0;
  for (iter = 0 ; iter < size ; iter++) {
    if (u_i[iter] == i && t_i[iter] == j) {
      exist = 1;
      break;
    }
  }
  return exist;
}

void GenWeights(int* unit_index, int* time_index, int* tr, int* C_it,
    int* len_data, /* total number of rows in the data */
    int* len_u_index,
    int* len_t_index,
    int* ate, int* att,
    double* weight) {
  int i;
  /* Rprintf("u_index: %d, t_index: %d, len_data: %d\n", *len_u_index, *len_t_index, *len_data); */
  
  for (i = 0 ; i < *len_u_index ; i++) {
    int j;
    for (j = 0 ; j < *len_t_index ; j++) {
      double w_it[*len_t_index];
      // initialize all elements in w_it to 0
      memset(w_it, 0, sizeof(double)*(*len_t_index));
      

      int c_it;
      int t_it;
      // initialize c_it and t_it
      int k;
      for (k = 0 ; k < *len_data ; k++) {
        if ( unit_index[k] == (i+1) && time_index[k] == (j+1) ) {
          c_it = C_it[k];
          t_it = tr[k];
          /* Rprintf(" t_it: %d\n", t_it); */
          /* Rprintf("0] unit: %d\n", i+1); */
        }
      }

  
      /* Rprintf("time index exist for unit %d and time %d: %d\n",
	 i+1, j+1, is_time_index_exist(unit_index, time_index, i+1,
	 j+1, *len_data )); */

      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
        if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
          // count control
          double count_control = 0;
          int k = 0;
          for (k = 0 ; k < *len_data ; k++) {
            if (unit_index[k] == (i+1) && tr[k] == 0)
              count_control++;
          }
          
          if (count_control > 0) { 
            double v_it = 1 / count_control;
            /* Rprintf("1] time: %d\n", j+1); */
            /* Rprintf("5.1] number of control: %f\n", count_control); */
            /* Rprintf("6] v_it: %f\n", v_it); */
            w_it[j] = 1;
            
            for (k = 0 ; k < *len_data ; k++) {
              if (unit_index[k] == (i+1) && tr[k] == 0) {
                int t_index = time_index[k] - 1;
                // Rprintf("2] opposite treatment index: %d\n", t_index+1);
                w_it[t_index] = v_it;
                
              }
            }

          }
          /* PdoubleArray(w_it, *len_t_index); */
          if (*ate == 1) {
            int k = 0;
            for (k = 0 ; k < *len_t_index ; k++)
	      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
		weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it);
	      }
          }
          else if (*att == 1) {
            int k = 0;
            for (k = 0 ; k < *len_t_index ; k++)
	      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
		weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it * t_it);
	      }
          }
        }        
        else if (t_it == 0) { /* ifTRUE(sub[,treat]...) == 0 */
          // count treate
          double count_treat = 0;
          int k = 0;
          for (k = 0 ; k < *len_data ; k++) {
            if (unit_index[k] == (i+1) && tr[k] == 1)
              count_treat++;
          }

          if (count_treat > 0) {
            double v_it = 1 / count_treat;
            /* Rprintf("3] j: %d\n", j+1); */
            /* Rprintf("5.2] number of treated: %f\n", count_treat); */
            /* Rprintf("6] v_it: %f\n", v_it); */
            w_it[j] = 1;
           
            for (k = 0 ; k < *len_data ; k++) {
              if (unit_index[k] == (i+1) && tr[k] == 1) {
                int t_index = time_index[k] - 1;
                /* Rprintf("4] opposite treatment index: %d\n", t_index+1); */
                w_it[t_index] = v_it;
              }
            }
          }
          /* PdoubleArray(w_it, *len_t_index); */
          if (*ate == 1) {
            int k = 0;
            for (k = 0 ; k < *len_t_index ; k++)
	      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
		weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it);
	      }
          }
          else if (*att == 1) {
            int k = 0;
            for (k = 0 ; k < *len_t_index ; k++)
	      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
		weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it * t_it);
	      }
          }

        }

      }
    }
  }
  
}


// Generate Weights for first-difference(FD)


void GenWeightsFD(int* unit_index, int* time_index, int* tr, int* C_it,
    int* len_data, /* total number of rows in the data */
    int* len_u_index,
    int* len_t_index,
    int* ate, int* att,
    double* weightfd) {
  int i;
  for (i = 0 ; i < *len_u_index ; i++) {
    /* Rprintf("1] unit index: %d\n", i); */
    int j;
    for (j = 0 ; j < *len_t_index ; j++) {
      double w_it[*len_t_index];
      // initialize all elements in w_it to 0
      memset(w_it, 0, sizeof(double)*(*len_t_index));
      
      int c_it;
      int t_it;
      // initialize t_it
      int k;
      for (k = 0 ; k < *len_data ; k++) {
        if ( unit_index[k] == (i+1) && time_index[k] == (j+1) ) {
	  c_it = C_it[k];          
	  t_it = tr[k];
          /* Rprintf("2] t_it: %d\n", t_it); */
        }
      }
      int t_it_1;
      // initialize t_it_1 (treatment for t-1)
      int m;
      for (m = 0 ; m < *len_data ; m++) {
        if ( unit_index[m] == (i+1) && time_index[m] == (j) ) {
          t_it_1 = tr[m];
          /* Rprintf("3] t_i(t-1): %d\n", t_it_1); */
        }
      }
      // checking whether previous time exit, i.e, j>1
      if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data) && is_time_index_exist(unit_index, time_index, i+1, j, *len_data)) {
        if ( t_it != t_it_1 ) { /* check whether treatment status changes  */
            w_it[j] = 1;
	    w_it[j-1] = 1;            
        }
        /* PdoubleArray(w_it, *len_t_index); */

	if (*ate == 1) {
	  int n = 0;
	  for (n = 0 ; n < *len_t_index ; n++)
	    if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
	      weightfd[n*(*len_u_index)+i] = weightfd[n*(*len_u_index)+i] + (w_it[n] * c_it);
	    }
	}
	else if (*att == 1) {
	  int n = 0;
	  for (n = 0 ; n < *len_t_index ; n++)
	    if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) {
	      weightfd[n*(*len_u_index)+i] = weightfd[n*(*len_u_index)+i] + (w_it[n] * c_it * t_it);
	    }
	}

      }
    }
  } 
}
