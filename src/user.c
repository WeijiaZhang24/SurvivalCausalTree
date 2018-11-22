/*
 * split.Rule = user (set temporarily as CT)
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"
#define SMLNUM (1e-3)

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;

int
userinit(int n, double *y[], int maxcat, char **error,
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



void
userss(int n, double *y[], double *value,  double *con_mean, double *tr_mean,
		double *risk, double *wt, double *treatment, double max_y,
		double alpha, double train_to_est_ratio, double *propensity, double *censoringProb, int *completeCase)
{
  //printf("\n%f max_y\n", max_y);
	int i;
	double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */
	double ttreat = 0.,tcon = 0.;
	double effect;
	double tr_var, con_var;
	double con_sqr_sum = 0., tr_sqr_sum = 0.;
	double pp = 0.;
	for (i = 0; i < n; i++) {
    //printf("%f\t",censoringProb[i]);
		temp1 += *y[i] * wt[i] * treatment[i] * completeCase[i] / propensity[i] / (censoringProb[i] + SMLNUM);
		temp0 += *y[i] * wt[i] * (1 - treatment[i]) * completeCase[i] / (1-propensity[i]) / (censoringProb[i] + SMLNUM);
    twt += wt[i];
		ttreat += wt[i] * treatment[i] / propensity[i] * completeCase[i] / (censoringProb[i] + SMLNUM);
		tcon += wt[i] * (1-treatment[i]) / (1-propensity[i]) * completeCase[i] / (censoringProb[i] + SMLNUM);
    //*/

    /*
		 temp1 += *y[i] * wt[i] * treatment[i] * completeCase[i] / propensity[i];
		 temp0 += *y[i] * wt[i] * (1 - treatment[i]) * completeCase[i] / (1-propensity[i]);
		 ttreat += wt[i] * treatment[i] / propensity[i] * completeCase[i] ;
		 tcon += wt[i] * (1-treatment[i]) / (1-propensity[i]) * completeCase[i] ;
		 //tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i] * completeCase[i] ;
		 //con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]) * completeCase[i];
		 //*/
    /*
		twt += wt[i];
		temp1 += *y[i] * wt[i] * treatment[i] / propensity[i];
		temp0 += *y[i] * wt[i] * (1 - treatment[i]) / (1-propensity[i]);
		ttreat += wt[i] * treatment[i] / propensity[i];
		tcon += wt[i] * (1-treatment[i]) / (1-propensity[i]);
		//tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		//con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
    //*/

	}

	effect = temp1 / ttreat - temp0 / tcon;
//	tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
//	con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));

	*tr_mean = temp1 / ttreat;
	*con_mean = temp0 / tcon;
	*value = effect;
	// this risk calculation doesn't make sense.
	//*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect;
	*risk = n-sqrt(n);

  // twt means "total weight"
	/*
	printf("\n%f\t value\n", *value);
	printf("%f\t risk\n", *risk);
	printf("%f\t tr_mean\n", *tr_mean);
	printf("%f\t con_mean\n", *con_mean);
	printf("%f\t twt\n", twt);
	printf("%f\t max_y\n", max_y);
	//*/
}


void user(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split,
		int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
		double train_to_est_ratio, double *propensity, double *censoringProb, int *completeCase)
{
	int i, j;
	double temp, temp2;
	double left_sum, right_sum;
	double left_tr_sum, right_tr_sum;
	double left_tr, right_tr;
	double left_wt, right_wt;
	int left_n, right_n;
	double best;
	int direction = LEFT;
	int where = 0;
	double node_effect, left_effect, right_effect;
	double left_temp, right_temp;
	int min_node_size = minsize;

	double tr_var, con_var;
	double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
	double left_tr_var, left_con_var, right_tr_var, right_con_var;

	//added definition
	double right_con_sum, right_con, right_tr2, left_con, left_con_sum, left_tr2;

	right_wt = 0.;
	right_tr = 0.;
	right_sum = 0.;
	right_tr_sum = 0.;
	right_sqr_sum = 0.;
	right_tr_sqr_sum = 0.;
	right_n = n;
	//added
	right_con_sum = 0.;
	right_con = 0.;
	right_tr2 = 0.;

	for (i = 0; i < n; i++) {
    //WJ: added

    right_tr2 += wt[i] * treatment[i]; // used for node size condition
		right_wt += wt[i]; // used for node size condition

    /*
		right_tr += wt[i] * treatment[i] / propensity[i];
		right_tr_sum += *y[i] * wt[i] * treatment[i] / propensity[i];
		right_con += wt[i] * (1-treatment[i]) / (1-propensity[i]);
		right_con_sum += *y[i] * wt[i] * (1-treatment[i]) / (1-propensity[i]);
		//*/

		///*
		right_tr += wt[i] * treatment[i] * completeCase[i] / propensity[i] / (censoringProb[i] + SMLNUM);
		right_con += wt[i] * (1-treatment[i]) * completeCase[i] / (1-propensity[i]) / (censoringProb[i] + SMLNUM);
		right_tr_sum += *y[i] * wt[i] * treatment[i] * completeCase[i]/ propensity[i] / (censoringProb[i]+ SMLNUM);
		right_con_sum += *y[i] * wt[i] * (1-treatment[i]) * completeCase[i] / (1-propensity[i]) / (censoringProb[i] +SMLNUM);
		//*/
	}

	temp = right_tr_sum / right_tr - right_con_sum / right_con;
	node_effect = alpha * temp * temp * right_wt;
  //printf("%f\n", node_effect);
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
		left_con_sum = 0;
		left_con = 0;
		left_tr2 = 0;
		for (i = 0; right_n > edge; i++) {
      //printf("%f\t", censoringProb[i]);
		  left_wt += wt[i];
		  right_wt -= wt[i];
		  left_tr2 += wt[i] * treatment[i];  // used for node size condition
		  right_tr2 -= wt[i] * treatment[i]; // used for node size condition

      left_tr += wt[i] * treatment[i] * completeCase[i] / propensity[i] / (censoringProb[i] + SMLNUM);
      left_con += wt[i] * (1-treatment[i]) * completeCase[i] / (1-propensity[i]) / (censoringProb[i] + SMLNUM);
      right_tr -= wt[i] * treatment[i] * completeCase[i] / propensity[i] / (censoringProb[i] + SMLNUM);
      right_con -= wt[i] * (1-treatment[i]) / (1-propensity[i]);

		  left_n++;
		  right_n--;

		  temp = *y[i] * wt[i] * treatment[i] * completeCase[i] / propensity[i] / (censoringProb[i] + SMLNUM);
		  left_tr_sum += temp;
		  right_tr_sum -= temp;
		  temp2 = *y[i] * wt[i] * (1-treatment[i]) * completeCase[i] / (1-propensity[i]) / (censoringProb[i] + SMLNUM);
		  left_con_sum += temp2;
		  right_con_sum -= temp2;
      /*
      printf("left_tr_sum\t %f\n", left_tr_sum);
      printf("left_tr\t %f\n", left_tr);
      printf("left_right_sum\t %f\n", right_tr_sum);
      printf("left_tr\t %f\n", right_tr);
      //*/
      //printf("above the loop\n");
		  if (x[i + 1] != x[i] && left_n >= edge &&
        (int) left_tr2 >= min_node_size &&
        (int) left_wt - (int) left_tr2 >= min_node_size &&
        (int) right_tr2 >= min_node_size &&
        (int) right_wt - (int) right_tr2 >= min_node_size) {
        //printf("in the loop\n");
		    left_temp = left_tr_sum / left_tr - left_con_sum / left_con;

		    right_temp = right_tr_sum / right_tr -right_con_sum / right_con;

		    left_effect = alpha * left_temp * left_temp * left_wt;
		    right_effect = alpha * right_temp * right_temp * right_wt;
		    temp = left_effect + right_effect - node_effect;
		    if (temp > best) {
		      best = temp;
		      where = i;
		      if (left_temp < right_temp)
		        direction = LEFT;
		      else
		        direction = RIGHT;
				}
			}
		}

		*improve = best;
		if (best > 0) {         /* found something */
			csplit[0] = direction;
			*split = (x[where] + x[where + 1]) / 2;
		}
	}
}


double
userpred(double *y, double wt, double treatment, double *yhat, double propensity)
{
	double ystar;
	double temp;

	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
	temp = ystar - *yhat;
	return temp * temp * wt;
}
