/*
 * split.Rule = CT
 */
#define OUTPUT 0
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"
static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;
int is_dominated(double *n1, double *n2, int k);

int
  USER4init(int n, double *y[], int maxcat, char **error,
            int *size, int who, double *wt, double *treatment,
            int bucketnum, int bucketMax, double *train_to_est_ratio)
  {
    if (who == 1 && maxcat > 0) {
      graycode_init0(maxcat);
      countn = (int *) ALLOC(2 * maxcat, sizeof(int));
      tsplit = countn + maxcat;
      treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
      wts = treatment_effect + maxcat;
      trs = wts + maxcat;
      sums = trs + maxcat;
      wtsums = sums + maxcat;
      trsums = wtsums + maxcat;
      wtsqrsums = trsums + maxcat;
      trsqrsums = wtsqrsums + maxcat;
    }
    *size = 1;
    *train_to_est_ratio = n * 1.0 / ct.NumHonest;
    return 0;
  }


void USER4ss(int n, double *y[], double *value,  double *con_mean, double *tr_mean,
             double *risk, double *wt, double *treatment, double max_y,
             double alpha, double train_to_est_ratio, double *propensity)
{
  int i;
  double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */
double ttreat = 0.;
double effect;
double tr_var, con_var;
double con_sqr_sum = 0., tr_sqr_sum = 0.;

for (i = 0; i < n; i++) {
  //temp1 += *y[i] * wt[i] * treatment[i] / propensity[i];
  temp1 += *y[i] * wt[i] * treatment[i];
  //temp0 += *y[i] * wt[i] * (1 - treatment[i]) / (1 - propensity[i]);
  temp0 += *y[i] * wt[i] * (1 - treatment[i]);
  twt += wt[i];
  ttreat += wt[i] * treatment[i];
  tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
  con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
}

effect = temp1 / ttreat - temp0 / (twt - ttreat);
tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));

*tr_mean = temp1 / ttreat;
*con_mean = temp0 / (twt - ttreat);
*value = effect;
*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect +
(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));
}


void USER4(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split,
           int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
           double train_to_est_ratio, double * propensity)
{
  int i, j;
  double temp, temp1, temp2, temp3;
  double left_sum, right_sum;
  double left_tr_sum, right_tr_sum;
  double left_tr, right_tr;
  double left_wt, right_wt;
  int left_n, right_n;
  double best, best1, best2, best3;
  int direction = LEFT;
  int where = 0, where1 = 0, where2 = 0, where3 = 0;
  int count;
  double node_effect, left_effect, right_effect;
  double left_temp, right_temp;
  int min_node_size = minsize;

  double tr_var, con_var;
  double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
  double left_tr_var, left_con_var, right_tr_var, right_con_var;

  //added definition
  double right_con_sum, right_con, right_tr2, left_con, left_con_sum, left_tr2;
  double node_effect2, left_effect2, right_effect2;


  right_wt = 0.;
  right_tr = 0.;
  right_sum = 0.;
  right_tr_sum = 0.;
  right_sqr_sum = 0.;
  right_tr_sqr_sum = 0.;
  right_n = n;



  //added
  double criterion1[right_n - edge], criterion2[right_n - edge], criterion3[right_n - edge];
  int split_index[right_n - edge];
  right_con_sum = 0.;
  right_con = 0.;
  right_tr2 = 0.;

  for (i = 0; i < n; i++) {
    right_wt += wt[i];
    right_tr2 += wt[i] * treatment[i]; // for node size condition
    right_tr += wt[i] * treatment[i];
    right_con += wt[i] * (1-treatment[i]);
    right_tr_sum += *y[i] * wt[i] * treatment[i];
    right_con_sum += *y[i] * wt[i] * (1-treatment[i]);

  }
  temp = right_tr_sum / right_tr - right_con_sum / right_con;
  node_effect = temp * temp * right_wt;
  node_effect2 = right_tr2 * (right_tr_sum / right_tr) * (right_tr_sum / right_tr) +
    right_con * (right_con_sum / right_con) * (right_con_sum / right_con);

  if (nclass == 0) {
    /* continuous predictor */
    left_wt = 0;
    left_tr = 0;
    left_n = 0;
    left_sum = 0;
    left_tr_sum = 0;
    left_sqr_sum = 0;
    left_tr_sqr_sum = 0;
    best = 0;
    //added
    best1 = 0;
    best2 = 0;
    best3 = 0;
    left_con_sum = 0;
    left_con = 0;
    left_tr2 = 0;
    count = 0;

    for (i = 0; right_n > edge; i++) {
      left_wt += wt[i];
      right_wt -= wt[i];
      left_tr2 += wt[i] * treatment[i];  // for node size condition
      right_tr2 -= wt[i] * treatment[i]; // for node size condition

      left_tr += wt[i] * treatment[i];
      left_con += wt[i] * (1-treatment[i]);
      right_tr -= wt[i] * treatment[i];
      right_con -= wt[i] * (1-treatment[i]);

      left_n++;
      right_n--;

      temp = *y[i] * wt[i] * treatment[i];
      left_tr_sum += temp;
      right_tr_sum -= temp;
      temp2 = *y[i] * wt[i] * (1-treatment[i]);
      left_con_sum += temp2;
      right_con_sum -= temp2;

      if (x[i + 1] != x[i] && left_n >= edge &&
          (int) left_tr2 >= min_node_size &&
          (int) left_wt - (int) left_tr2 >= min_node_size &&
          (int) right_tr2 >= min_node_size &&
          (int) right_wt - (int) right_tr2 >= min_node_size) {

        left_temp = left_tr_sum / left_tr - left_con_sum / left_con;
        right_temp = right_tr_sum / right_tr -right_con_sum / right_con;

        left_effect = left_temp * left_temp * left_wt;
        right_effect = right_temp * right_temp * right_wt;

        left_effect2 = left_tr2 * (left_tr_sum / left_tr) * (left_tr_sum / left_tr)
          + (left_wt - left_tr2) * (left_con_sum / left_con) * (left_con_sum / left_con);
        right_effect2 = right_tr2 * (right_tr_sum / right_tr) * (right_tr_sum / right_tr)
          + right_con * (right_con_sum / right_con) * (right_con_sum / right_con);

        //temp1 = fabs(left_temp - right_temp);			  		//MH
        temp1 = left_temp * left_temp + right_temp * right_temp + 2 * left_temp * right_temp;

        temp3 = left_effect + right_effect - node_effect; 		//CT
        temp2 = left_effect2 + right_effect2 - node_effect2; 	//fit

        criterion1[count] = temp1;
        criterion2[count] = temp2;
        criterion3[count] = temp3;
        split_index[count] = i;
        //track how many times go inside this loop
        count++;

        if (temp1 > best1) {
          best1 = temp1;
          where1 = i;
          if (left_temp < right_temp)
            direction = LEFT;
          else
            direction = RIGHT;
        }
        if (temp2 > best2) {
          best2 = temp2;
          where2 = i;
          if (left_temp < right_temp)
            direction = LEFT;
          else
            direction = RIGHT;
        }
        if (temp3 > best3) {
          best3 = temp3;
          where3 = i;
          if (left_temp < right_temp)
            direction = LEFT;
          else
            direction = RIGHT;
        }
        //*/
      }
    }

    int num_dominated[count];

    double a[3], b[3];
    for (i=0; i < count; i++)
    {
      num_dominated[i]=0;
      for (j=0; j < count; j++)
      {
        a[0] = criterion1[i]; a[1] = criterion2[i]; a[2] = criterion3[i];
        b[0] = criterion1[j]; b[1] = criterion2[j]; b[2] = criterion3[j];
        if (is_dominated(a,b,2) ==1 || criterion2[i] <= 0 || criterion1[i] <=0 || criterion3[i] <= 0)
          num_dominated[i]++;
      }
    }
    if (where1 != 0 && OUTPUT){
      printf("MHindex \t %d \t MH \t %f \n",where1, best1);
      printf("CTindex \t %d \t CT \t %f \n", where2, best2);
      printf("FITindex \t %d \t FIT \t %f \t \n", where3, best3);
    }
    if (count > 0){
      for (i =0; i < count; i++){
        if (num_dominated[i]==0)
        {
          // first try random
          /*
          int n = rand() % 2;
          if (where == 0 || n ==1){
          where = split_index[i];
          best = criterion1[i];
          }
          // */
          // no strategy, just use the last one
            where = split_index[i];
            best = criterion2[i];
        }
      }
  }
    else{
      best = 0;
      where = 0;
    }
    if (where != 0 && OUTPUT)
      printf("out of loop\n");
    direction = LEFT;
    *improve = best;
    if (best > 0) {         /* found something */
    csplit[0] = direction;
      *split = (x[where] + x[where + 1]) / 2;
    }
}
  }




double
  USER4pred(double *y, double wt, double treatment, double *yhat, double propensity)
  {
    double ystar;
    double temp;

    ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
    temp = ystar - *yhat;
    return temp * temp * wt;
  }

