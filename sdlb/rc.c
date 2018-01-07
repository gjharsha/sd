//
//  rc.c
//  sdlb
//
//  Created by Yifan Liu on 6/26/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//


#include <time.h>
#include <float.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "theta.h"
#include "optimal.h"
#include "omega.h"
#include "log.h"
#include "master.h"
#include "cuts.h"
#include "rvgen.h"
#include "sdconstants.h"
#include "sdglobal.h"
#include "resumeb.h"
#include "lambda.h"
#include "sigma.h"
#include "delta.h"
#include "rc.h"

BOOL get_index_number(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s, int obs_idx)
{
    int i,j;
    int *cstat;
    int *rstat;
    char *basismsg;
    unsigned long *plain, *plain2;
    BOOL New_index=TRUE;
    
    /* One extra bit for accomodating index and another for one norm */
    /* Save the 0th location for the norm of the index vector */
    if (!(plain = arr_alloc(s->ids->num_word, unsigned long)))
		err_msg("Allocation", "get_index_number", "plain");
    /* plain2 is for basis status of the slack */
    if (!(plain2 = arr_alloc(s->ids->num_word2, unsigned long)))
        err_msg("Allocation", "get_index_number", "plain2");
    
	cstat = (int *) malloc((p->num->sub_cols + 1) * sizeof(int));
    rstat = (int *) malloc((p->num->sub_rows + 1) * sizeof(int));
    
    /* First of all, let record the current obeservation idx */
    s->ids->current_obs_idx = obs_idx;
    
	get_basis(c->subprob, cstat + 1, rstat + 1); /* 2011.10.30 */
    
	if (0)
	{
		/* Write out the solution */
        double *y;
        y = (double *) malloc((p->num->sub_cols + 1) * sizeof(double));
        get_x(c->subprob, y + 1, 0, p->num->sub_cols - 1); /* 2011.10.30 */
		for (j = 1; j <= p->num->sub_cols; j++)
		{
			printf("y[%d] = %17.10g", j, y[j]);
			if (cstat != NULL)
			{
				switch (cstat[j])
				{
                    case AT_LOWER:
                        basismsg = "Nonbasic at lower bound";
                        break;
                    case BASIC:
                        basismsg = "Basic";
                        break;
                    case AT_UPPER:
                        basismsg = "Nonbasic at upper bound";
                        break;
                    case FREE_SUPER:
                        basismsg = "Superbasic, or free variable at zero";
                        break;
                    default:
                        basismsg = "Bad basis status";
                        break;
				}
				printf("  %s   %d", basismsg, cstat[j]);
			}
			printf("\n");
		}
        printf("obj=:%f\n",get_objective(c->subprob));
        write_prob(c->subprob, "deg_sub_check.lp");
	}
    
    if (0) {
        for (j = 1; j <= p->num->sub_cols; j++) {
            printf("%d", cstat[j]);
        }
        printf("\t %d\n",c->k);
    }

    
#if CHECK_PHI
    int countc=0,countr=0;
    for (j = 1; j <= p->num->sub_cols; j++) {
        countc = countc + cstat[j];
    }
    printf("countc = %d\t",countc);
    for (j = 1; j <= p->num->sub_rows; j++) {
        countr = countr + rstat[j];
    }
    printf("countr = %d\t sum=%d\n",countr,countc+countr);
#endif

    encode_col(p, plain, cstat, WORD_LENGTH);
    for (j = 1; j < s->ids->num_word; j++) {
        plain[0] += plain[j];
    }
    encode_row(p, plain2, rstat, WORD_LENGTH);
    for (j = 1; j < s->ids->num_word2; j++) {
        plain2[0] += plain2[j];
    }
    
    for (i = 0; i < s->ids->cnt; i++) {
        if (equal_ulong_arr(plain, s->ids->index[i]->val, s->ids->num_word) && equal_ulong_arr(plain2, s->ids->index2[i]->val, s->ids->num_word2)) {
            New_index = FALSE;
            s->ids->index[i]->freq++;
            s->ids->current_index_idx = i;
        }
    }
    
    if (New_index) {
        s->ids->index[s->ids->cnt] = new_id(s->ids->num_word);
        s->ids->index[s->ids->cnt]->val = plain;
        s->ids->index[s->ids->cnt]->first_c_k = c->k;
        s->ids->index[s->ids->cnt]->freq = 1;
        s->ids->index2[s->ids->cnt] = new_id(s->ids->num_word2);
        s->ids->index2[s->ids->cnt]->val = plain2;
        s->ids->index2[s->ids->cnt]->first_c_k = c->k;
        s->ids->index2[s->ids->cnt]->freq = 1;
        s->ids->current_index_idx = s->ids->cnt;
        s->ids->cnt++;
    }
    
#ifdef CHECK_PHI
    printf("%lu\n",s->ids->omega_index[1]);
#endif
    
    if (p->num->rv_g && New_index) {
        
        /* modified by Yifan 2014.08.01 There's no need to do bitwise operation if there is no random cost */
        for (j = 1; j<s->ids->num_word; j++) {
            /* Bitwise "B and O" */
            s->rcdata->phi_col[j] = plain[j] & s->ids->omega_index[j];
            /* Bitwise "(not B) and O " */
            s->rcdata->lhs_chl[j] = (~plain[j]) & s->ids->omega_index[j];
        }
        
        /* modified by Yifan 2014.08.01 There's no need to decode anything if there is no random cost*/
#ifdef CHECK_PHI
        printf("Here are the Basis columns:\n");
#endif
        /* clean up col_num first */
        for (j = 0; j < p->num->sub_cols; j++) {
            s->rcdata->col_num[j] = 0;
        }
        decode_col(p, s->rcdata->col_num, plain, WORD_LENGTH);
        /* Let's decode Phi column. Note: phi_col's zero-th location store the actual norm */
#ifdef CHECK_PHI
        printf("Here are the Phi columns:\n");
#endif
        /* clean up phi_col_num first */
        for (j = 0; j < p->num->rv_g; j++) {
            s->rcdata->phi_col_num[j] = 0;
        }
        s->rcdata->phi_col[0] = decode_col(p, s->rcdata->phi_col_num, s->rcdata->phi_col, WORD_LENGTH);
        if (s->rcdata->phi_col[0]) {
            new_phi(s,p,(int) s->rcdata->phi_col[0]);
            get_phi_val(s, p, s->ids->index[s->ids->current_index_idx]);
            get_cost_val(sd_global, s->omega, p->num, s->ids->index[s->ids->current_index_idx], s->ids->random_cost_val, s->ids->random_cost_col, s->ids->current_obs_idx);
        }
        else{
            s->ids->index[s->ids->current_index_idx]->phi_cnt = (int) s->rcdata->phi_col[0];
        }
        /* Let's decode lhs column. Note: lhs_col's zero-th location store the actual norm */
#ifdef CHECK_PHI
        printf("Here are the Nonbasis random columns:\n");
#endif
        s->rcdata->lhs_chl[0] = decode_col(p, s->rcdata->lhs_col_num, s->rcdata->lhs_chl, WORD_LENGTH);

    }

#if 0
    FILE *index_number;
    index_number = fopen("indexNumber.txt", "a");
    for (j = 1; j <= p->num->sub_cols; j++) {
            fprintf(index_number, "%d\t", cstat[j]);
        }
    fprintf(index_number, "\n");
    fclose(index_number);
#endif
    
    if (!New_index) {
        mem_free(plain);
        mem_free(plain2);
    }
    
    mem_free(cstat);
    mem_free(rstat);
    
    return New_index;
}
/* Decode the column and return the total number of 1's */
int decode_col(prob_type *p, int *col_num, unsigned long *col, int word_length)
{
    /* Use unsigned long, so that when shift the bit right, left will be padded with 0's */
    unsigned long mask, temp;
    mask = 1;
    int j, cnt, group, shift;
    /* Let's decode phi_col */
    cnt = 0;
    for (j = 1; j <= p->num->sub_cols; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) col[group] >> shift;
        temp = temp & mask;
        if (temp) {
            col_num[cnt] = j;
#ifdef CHECK_PHI
            printf("column number is %d\n", col_num[cnt]);
#endif
            cnt++;
        }
    }
    if (!cnt) {
#ifdef CHECK_PHI
        printf("None!\n");
#endif
    }
    
    return cnt;
}

int encode_col(prob_type *p, unsigned long *col, int *cstat, int word_length)
{
    int j,group,shift;
    unsigned long temp;
    for (j = 1; j <= p->num->sub_cols; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) cstat[j] << shift;
        col[group] |= temp;
    }
    return 0;
}

int encode_row(prob_type *p, unsigned long *row, int *rstat, int word_length)
{
    int j,group,shift;
    unsigned long temp;
    for (j = 1; j <= p->num->sub_rows; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) rstat[j] << shift;
        row[group] |= temp;
    }
    return 0;
}

/***********************************************************************\
 ** This function loops through all the dual vectors found so far
 ** and returns the index of the one which satisfies the expression:
 **
 **	argmax { Pi x (R - T x X) | all Pi }
 **
 ** where X, R, and T are given.  It is calculated in this form:
 **
 ** 	Pi x Rbar + Pi x Romega + (Pi x Tbar) x X + (Pi x Tomega) x X
 **
 ** Since the Pi's are indexed by the index sets which point to two
 ** different structures (sigma and delta). The following procedure goes
 ** through all the index set and retrieves all sigma and delta needed
 ** for the calculation.
 \***********************************************************************/
i_type compute_istar_index(sdglobal_type *sd_global, soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                     delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                     double *argmax, BOOL pi_eval, int ictr)
{
    double arg;
    //  double	argmax;             //modified by Yifan to return argmax value 09/22/2011
    int sig_pi, del_pi, index_idx;
    int c, new_pisz;
    i_type ans;
    
#ifdef LOOP
    printf("Inside compute_istar_index\n");
#endif
    //printf("\n***IN COMPUTE ISTAR ...\n");
    if (pi_eval == TRUE)
    {
        //new_pisz =sd_global->config.PI_EVAL_START;
        new_pisz = ictr / 10 + 1;
    }
    else
        new_pisz = 0;
    
    ictr -= new_pisz; /*evaluate the pi's generated in the first 90% iterations */
    
    *argmax = -DBL_MAX;
    
    /*added by Yifan to enable parallel process*/
    // #pragma omp for private(sig_pi,del_pi, arg)
    for (index_idx = 0 ; index_idx < s->ids->cnt; index_idx++) {
        sig_pi = s->ids->sig_idx[index_idx];
        if (sigma->ck[sig_pi] <= ictr) {
            /* Find the row in delta corresponding to this row in sigma */
            del_pi = sigma->lamb[sig_pi];
            
            /* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
            arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R
            - Pi_Tbar_X[sig_pi];
            
            /* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
            
            for (c = 1; c <= num->rv_cols; c++)
                arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];
            
            if (num->rv_g && s->ids->index[index_idx]->phi_cnt) {
                get_cost_val(sd_global, s->omega, num, s->ids->index[index_idx], s->ids->random_cost_val, s->ids->random_cost_col, obs);
                adjust_argmax_value(s, obs, sigma, delta, Xvect, num, Pi_Tbar_X, &arg, s->ids->index[index_idx]);
            }
            
#ifdef LOOP
            print_sigma(sigma, num, sig_pi);
            print_delta(delta, num, del_pi, obs);
            printf("\nResulting arg=%f", arg);
#endif
            
            if (arg > (*argmax))
            {
                *argmax = arg;
                ans.sigma = sig_pi;
                ans.delta = del_pi;
                ans.index_idx = index_idx;
                //printf("argmax:%f and istar(%d,%d)\n", arg, ans.sigma, ans.delta);
            }
        }
    }
    
#ifdef LOOP
    printf("Exiting compute_istar_index\n");
#endif
    
    return ans;
}

i_type compute_new_istar_index(sdglobal_type *sd_global, soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                         delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                         double *argmax, int ictr)
{
    double arg;
    int sig_pi, del_pi, index_idx;
    int c, new_pisz;
    i_type ans;
    ans.sigma = 0;
    ans.delta = 0;
    
#ifdef LOOP
    printf("Inside compute_new_istar_index\n");
#endif
    
    new_pisz = ictr / 10 + 1;
    ictr -= new_pisz; /*evaluate the pi's generated in the last 10% iterations */
    
    *argmax = -DBL_MAX;
    
    for (index_idx = 0; index_idx < s->ids->cnt; index_idx++) {
        sig_pi = s->ids->sig_idx[index_idx];
        if (sigma->ck[sig_pi] > ictr) {
            /* Find the row in delta corresponding to this row in sigma */
            del_pi = sigma->lamb[sig_pi];
            
            /* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
            arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R - Pi_Tbar_X[sig_pi];
            
            /* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
            for (c = 1; c <= num->rv_cols; c++)
                arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];
            
            if (num->rv_g && s->ids->index[index_idx]->phi_cnt) {
                get_cost_val(sd_global, s->omega, num, s->ids->index[index_idx], s->ids->random_cost_val, s->ids->random_cost_col, obs);
                adjust_argmax_value(s, obs, sigma, delta, Xvect, num, Pi_Tbar_X, &arg, s->ids->index[index_idx]);
            }
#ifdef LOOP
            print_sigma(sigma, num, sig_pi);
            print_delta(delta, num, del_pi, obs);
            printf("\nResulting arg=%f", arg);
#endif
            
            if (arg > (*argmax))
            {
                *argmax = arg;
                ans.sigma = sig_pi;
                ans.delta = del_pi;
                ans.index_idx = index_idx;
            }
        }
    }
    
#ifdef LOOP
    printf("Exiting compute_new_istar_index\n");
#endif
    return ans;
}

void new_phi(soln_type *s, prob_type *p, int num_of_phi)
{
    int i;
    if (!(s->ids->index[s->ids->current_index_idx]->phi_val = mem_calloc(num_of_phi, sizeof(double *))))
        err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_val");
    for (i = 0; i < num_of_phi; i++) {
        if (!(s->ids->index[s->ids->current_index_idx]->phi_val[i] = arr_alloc(p->num->sub_rows+1, double)))
            err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_val");
    }
    if (!(s->ids->index[s->ids->current_index_idx]->phi_cost_delta = arr_alloc(num_of_phi, double)))
        err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_cost_delta");
    if (!(s->ids->index[s->ids->current_index_idx]->phi_col_num = arr_alloc(num_of_phi, int)))
        err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_val");
    if (!(s->ids->index[s->ids->current_index_idx]->phi_sigma_idx = arr_alloc(num_of_phi, int)))
        err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_sigma_idx");
    if (!(s->ids->index[s->ids->current_index_idx]->phi_lambda_idx = arr_alloc(num_of_phi, int)))
        err_msg("Allocation", "new_phi", "s->ids->index[s->ids->current_index_idx]->phi_lambda_idx");
    
    s->ids->index[s->ids->current_index_idx]->phi_cnt = num_of_phi;

#if 0
    printf("\tcurrent_idx:%d;\tnum_of_phi:%d\n",s->ids->current_index_idx, num_of_phi);
#endif
}

void free_phi(id_type *index)
{
    int i;
    for (i = 0; i < index->phi_cnt; i++) {
        mem_free(index->phi_val[i]);
    }
    mem_free(index->phi_val);
    mem_free(index->phi_cost_delta);
    mem_free(index->phi_col_num);
    mem_free(index->phi_sigma_idx);
    mem_free(index->phi_lambda_idx);
}

int get_phi_val(soln_type *s, prob_type *p, id_type *index)
{
    int *newhead;
    double *phi;
    int i,j,k,idx=0;
    if (!(newhead = mem_calloc(p->num->sub_rows, sizeof(int))))
        err_msg("Allocation", "get_phi_val", "newhead");
    if (!(phi = mem_calloc(p->num->sub_rows, sizeof(double))))
        err_msg("Allocation", "get_phi_val", "phi");
    
    /* Let's print out all the phi value */
    
    /* First, get the mapping of basis from solver */
    get_basis_head(p->subprob, newhead);
    
#if 0
    /* Check the content of the basis mapping */
    for (i = 0; i < p->num->sub_rows; i++) {
        printf("newhead[%d]=%d\n", i, newhead[i]);
    }
#endif
    
    for (i = 0; i < p->num->sub_rows; i++) {
        for (j = 0; j < p->num->rv_g; j++) {
            
            /* Note the newhead starts with 0 while s->rcdata->phi_col_num starts with 1 */
            if (newhead[i] == s->rcdata->phi_col_num[j]-1) {
                /* Get phi and save a copy */
                get_basis_row(p->subprob, i, phi);
                copy_arr(index->phi_val[idx]+1, phi, p->num->sub_rows-1);
                /* Calculate the one norm of phi */
                index->phi_val[idx][0] = one_norm(index->phi_val[idx]+1, p->num->sub_rows);
                /* Save the column number of phi */
                index->phi_col_num[idx] = s->rcdata->phi_col_num[j];
                
#if 0
                /* Check the content of phi */
                printf("The following is the phi we want for column %d:\n", index->phi_col_num[idx]);
                for(k = 0; k < p->num->sub_rows; k++) {
                    printf("phi[%d]:%f\n", k+1, phi[k]);
                }
                printf("\n");
#endif
                idx++;
            }
        }
    }
    
    mem_free(newhead);
    mem_free(phi);
    return 0;
}

/* This function obtains the value of each cost coefficient and the corresponding column number */
int get_cost_val(sdglobal_type *sd_global, omega_type *omega, num_type *num, id_type *index, double *cost, int *col, int obs_idx)
{
    int cnt, idx;
    sparse_vect Gomega;
    init_G_omega(&Gomega, omega, num);
    get_G_omega(sd_global, num, omega, obs_idx);
    for (cnt = 1; cnt <= Gomega.cnt; cnt++)
    {
        cost[cnt] = Gomega.val[cnt];
        col[cnt] = Gomega.row[cnt]-num->mast_cols;
        for (idx = 0; idx < index->phi_cnt; idx++) {
            if (col[cnt] == index->phi_col_num[idx]) {
                /* phi_cost_delta starts its storage from index zero. It will match the phi_val's starting location */
                index->phi_cost_delta[idx] = cost[cnt];
            }
        }
    }

    return 0;
}
/* If random cost exists and it will affect the dual solution, then we need to calculate Nu.
   Pi = Nu + sum_{j} cost_delta(j)*phi(j)
  delta_g(j) is a scalar and theta(j) is a column vector */
int adjust_dual_solution(double *Pi, num_type *num, id_type *index)
{
    int i,j;
    
    if (index->phi_cnt) {
        for (i = 1; i <= num->sub_rows; i++) {
            for (j = 0; j < index->phi_cnt; j++) {
                Pi[i] -= index->phi_cost_delta[j] * index->phi_val[j][i];
            }
        }
    }
    Pi[0] = one_norm(Pi + 1, num->sub_rows);
    return 0;
}

int put_phi_into_sigma_delta(sdglobal_type* sd_global, cell_type *c, soln_type *s,
                             lambda_type *lambda, sigma_type *sigma, delta_type *delta,
                             omega_type *omega, num_type *num, sparse_vect *Rbar,
                             sparse_matrix *Tbar,id_type *index)
{
    int cnt;
    int lamb_idx;
    BOOL new_lamb = FALSE, new_sigma = FALSE;
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        lamb_idx = calc_lambda(sd_global, lambda, num, index->phi_val[cnt], &new_lamb);
        index->phi_lambda_idx[cnt] = lamb_idx;
        index->phi_sigma_idx[cnt] = calc_sigma(sd_global, c, sigma, num, index->phi_val[cnt], Rbar, Tbar, lamb_idx, new_lamb, &new_sigma);
        if (new_lamb) {
            calc_delta_row(sd_global, delta, lambda, omega, num, lamb_idx);
        }
    }
    return 0;
}

int adjust_argmax_value(soln_type *s, int obs, sigma_type *sigma, delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,double *arg, id_type *index)
{
    int cnt;
    int sig_pi, del_pi, c;
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        sig_pi = index->phi_sigma_idx[cnt];
        del_pi = sigma->lamb[sig_pi];
        
        /* First need to update the phi_cost_delta value according obs */
        
        /* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
        *arg += index->phi_cost_delta[cnt] * (sigma->val[sig_pi].R + delta->val[del_pi][obs].R - Pi_Tbar_X[sig_pi]);
        /* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
        for (c = 1; c <= num->rv_cols; c++) {
            *arg -= index->phi_cost_delta[cnt] * delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];
        }
    }
    return 0;
}

int adjust_alpha_value(soln_type *s, int obs, one_cut *cut ,sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index, int weight)
{
    int cnt;
    int sig_pi, del_pi;
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        sig_pi = index->phi_sigma_idx[cnt];
        del_pi = sigma->lamb[sig_pi];
        cut->alpha += index->phi_cost_delta[cnt] * sigma->val[sig_pi].R * weight;
        cut->alpha += index->phi_cost_delta[cnt] * delta->val[del_pi][obs].R * weight;
    }
    
    return 0;
}

int adjust_beta_value(soln_type *s, int obs, one_cut *cut ,sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index, int weight)
{
    int cnt;
    int sig_pi, del_pi, c;
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        sig_pi = index->phi_sigma_idx[cnt];
        del_pi = sigma->lamb[sig_pi];
        for (c = 1; c <= num->nz_cols; c++)
            cut->beta[sigma->col[c]] += index->phi_cost_delta[cnt] * sigma->val[sig_pi].T[c]
            * weight;
        for (c = 1; c <= num->rv_cols; c++)
            cut->beta[delta->col[c]] += index->phi_cost_delta[cnt] * delta->val[del_pi][obs].T[c]
            * weight;
    }
    return 0;
}

int adjust_incumbent_height(soln_type *s, int obs, one_cut *cut , double *beta, sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index)
{
    int cnt, c;
    
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        cut->subobj_omega[obs] += index->phi_cost_delta[cnt] * sigma->val[index->phi_sigma_idx[cnt]].R;
        cut->subobj_omega[obs] += index->phi_cost_delta[cnt] * delta->val[index->phi_lambda_idx[cnt]][obs].R;
    }
    
    for (cnt = 0; cnt < index->phi_cnt; cnt++) {
        for (c = 1; c <= num->nz_cols; c++)
            beta[sigma->col[c]] += index->phi_cost_delta[cnt] * sigma->val[index->phi_sigma_idx[cnt]].T[c];
        for (c = 1; c <= num->rv_cols; c++)
            beta[delta->col[c]] += index->phi_cost_delta[cnt] * delta->val[index->phi_lambda_idx[cnt]][obs].T[c];
    }
    return 0;
}
