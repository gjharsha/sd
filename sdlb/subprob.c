/***********************************************************************\
**
 ** subprob.c
 **
 ** This file contains the functions needed for solving the subproblem
 **
 **
 **
 ** solve_subprob()
 ** compute_rhs_cost()
 **
 ** History:
 **   20 Jan 1992 - <Jason Mai> - created.
 **   ?? Feb 1992 - <Jason Mai> - updated with new structures.
 **   22 Mar 1992 - <Jason Mai> - debugged. 
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include <time.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "testout.h"
#include "subprob.h"
#include "omega.h"
#include "log.h"
#include "solver.h"
#include "sdglobal.h"

/***********************************************************************\
** This function will solve a new subproblem. This involves
 ** replacing the right hand side of the subproblem with new values,
 ** based upon some observation of omega, and some X vector of primal
 ** variables from the master problem.  Generally, the latest observation
 ** is used.  When forming a normal cut, the candidate x should be used, while
 ** the incumbent x should be used for updating the incumbent cut.
 \***********************************************************************/
/* */
int solve_subprob(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, vector Xvect, int omeg_idx)
{
	vector rhs;
    vector cost;
	BOOL ans;
	clock_t start, end; /* Recording solution time for solving
	 subproblem LPs. added by zl, 06/29/04. */
	double sub_obj; /* To check the valid LB on subproblem
	 objective function values. zl, 07/01/04. */

#ifdef SAVE
	char fname[20] = "subprob   .lp";
	static int fnum = 0;
#endif

#ifdef TRACE
	printf("Inside solve_subprob\n");
#endif
    if (!(rhs = arr_alloc(p->num->sub_rows+1, double)))
		err_msg("Allocation", "compute_rhs", "rhs");
    if (!(cost = arr_alloc(p->num->sub_cols+1, double)))
		err_msg("Allocation", "compute_rhs", "rhs");

	compute_rhs_cost(sd_global, p->num, p->Rbar, p->Tbar, p->Gbar, Xvect, rhs, cost, s->omega, omeg_idx);
	if (!change_col(c->subprob, RHS_COL, rhs + 1, 0, p->num->sub_rows))
	{
		print_contents(c->subprob, "contents.out");
		err_msg("change_col", "solve_subprob", "returned FALSE");
	}
    
    /* modified by Yifan 2014.06.11 */
    if (p->num->rv_g > 0) {
        if (!change_row(c->subprob, -1, cost + 1, 0, p->num->sub_cols))
        {
            print_contents(c->subprob, "contents.out");
            err_msg("change_col", "solve_subprob", "returned FALSE");
        }
    }
    /* print_problem(c->subprob, "testing_sub.lp"); */
    
    /* modified by Yifan 2014.06.17 */
#if 0
    FILE *rhs_value;
    rhs_value =fopen("rhs_value.txt", "a");
    fprintf(rhs_value, "%f\t%f\n",rhs[1],rhs[2]);
    fclose(rhs_value);
    FILE *cost_value;
    cost_value =fopen("cost_value.txt", "a");
    fprintf(cost_value, "%f\n",cost[3]);
    fclose(cost_value);
#endif

#ifdef DEBUG
	print_vect(rhs, p->num->sub_rows, "Rhs");
#endif

#ifdef SAVE
	fname[7] = '0' + fnum / 100 % 10;
	fname[8] = '0' + fnum / 10 % 10;
	fname[9] = '0' + fnum / 1 % 10;
	++fnum;
	print_problem(c->subprob, fname);
	printf("Saving file: %s\n", fname);
#endif


	/* Recording the time for solving subproblem LPs. zl, 06/29/04. */
	start = clock();
	c->subprob->feaflag = TRUE; /*added by Yifan to generate feasibility cut 08/11/2011*/
	ans = solve_problem(sd_global, c->subprob);
	end = clock();
	s->run_time->soln_subprob_iter += ((double) (end - start)) / CLOCKS_PER_SEC;

	/* Record the lowest subproblem objective function values so far.
	 zl, 07/01/04. */
	sub_obj = get_objective(c->subprob);

	if (s->sub_lb_checker > sub_obj)
		s->sub_lb_checker = sub_obj;

	c->LP_cnt++; /* # of LPs solved increase by 1. zl 06/30/02 */
    
    mem_free(rhs);
    mem_free(cost);

#ifdef TRACE
	printf("Exiting solve_subprob\n");
#endif

	return ans;
}

/***********************************************************************\
** This function computes the right hand side of the subproblem,
 ** based upon a given X vector and a given observation of omega.
 ** It is defined as:
 **	rhs = R(omega) - T(omega) x X
 **
 ** and is calculated as:
 ** 	rhs = Rbar - Tbar x X + Romega - Tomega x X
 **
 ** where the "bar" denotes the fixed or mean value, and the "omega"
 ** denotes a random variation from this mean.
 ** The function allocates an array for the vector, which must be
 ** freed by the customer.  Also, the zeroth position of this rhs vector
 ** is reserved, and the actual values begin at rhs[1].
 \***********************************************************************/
void compute_rhs_cost(sdglobal_type* sd_global, num_type *num, sparse_vect *Rbar,
		sparse_matrix *Tbar, sparse_vect *Gbar, vector X, vector rhs, vector cost, omega_type *omega, int omeg_idx)
{
	int cnt;
	sparse_vect Romega;
	sparse_matrix Tomega;
    sparse_vect Gomega;

#ifdef TRACE
	printf("Inside compute_rhs_cost\n");
#endif

	init_R_T_G_omega(&Romega, &Tomega, &Gomega, omega, num);
	get_R_T_G_omega(sd_global, omega, omeg_idx);

#ifdef DEBUG
	printf("\n\n");
	print_vect(X, num->mast_cols, "Candidate X");
	print_sparse_vect(Rbar, "Rbar");
	print_sparse_vect(&Romega, "Romega");
	print_sparse_matrix(Tbar, "Tbar");
	print_sparse_matrix(&Tomega, "Tomega");
	printf("\n\n");
#endif

	for (cnt = 0; cnt <= num->sub_rows; cnt++)
		rhs[cnt] = 0.0;

	/* Start with the values of R(omega) -- both fixed and varying */
	for (cnt = 1; cnt <= Rbar->cnt; cnt++)
		rhs[Rbar->row[cnt]] += Rbar->val[cnt];
	for (cnt = 1; cnt <= Romega.cnt; cnt++)
		rhs[Romega.row[cnt]] += Romega.val[cnt];

	/* (cumulatively) subtract values of T(omega) x X -- both fixed and varying */
	TxX(Tbar, X, rhs);
	TxX(&Tomega, X, rhs);

	/* This rhs vector is one element too large. */
    
    /* modified by Yifan 2014.06.11 */
    /* End with the values of G(omega) -- both fixed and varying */
	for (cnt = 1; cnt <= Gbar->cnt; cnt++)
		cost[Gbar->row[cnt]] += Gbar->val[cnt];
	for (cnt = 1; cnt <= Gomega.cnt; cnt++)
		cost[Gomega.row[cnt]-num->mast_cols] += Gomega.val[cnt];

#ifdef TRACE
	printf("Exiting compute_rhs_cost\n");
#endif
}

/***********************************************************************\
** Once every cell needs its own copy of the subproblem, this function
 ** will be called from solve_cell to make a duplicate of the pristine
 ** subproblem stored in the prob structure.  Until then, every cell
 ** can just re-use the same subproblem, and this function just returns
 ** the pointer it was passed.
 \***********************************************************************/
one_problem *new_subprob(one_problem *subprob)
{
#ifdef TRACE
	printf("Inside new_subprob\n");
#endif

	if (!(setup_problem(subprob)))
		err_msg("Problem Setup", "new_subprob", "subprob");

	return subprob;
}

/***********************************************************************\
** Since, for the serial case and the simulated parallel case, the 
 ** subproblem is not actually copied, we don't want to free it here.  
 ** We just remove it from CPLEX so the next cell can load it back in.
 \***********************************************************************/
void free_subprob(one_problem *subprob)
{
#ifdef TRACE
	printf("Inside free_subprob\n");
#endif

	remove_problem(subprob);
}

