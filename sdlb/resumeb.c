//
//  resumeb.c
//  sdlb
//
//  Created by Yifan Liu on 1/12/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//

#include <time.h>
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

int store_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int idx, obs, cnt;
    char buffer1[128], buffer2[128];
    FILE *rep_data;
    double coef_value;
    
    char rep_number[16];
    sprintf(rep_number, "%d", p->current_batch_id);
    
    /* modified by Yifan 2014.02.26 Numbering resume data for different replications */
    strcpy(buffer2, "resume");
    strcat(buffer2, rep_number);
    strcat(buffer2, ".lp");
    write_prob(c->master,buffer2);
    
    
    /* modified by Yifan 2014.02.26 */
    strcpy(buffer1, "resume_data");
    strcat(buffer1, rep_number);
    strcat(buffer1, ".txt");
    
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen(buffer1, "w");
    
    
    /* 1. Start storing omega_type */
    fwrite(&(s->omega->cnt), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->next), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->most), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->last), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->k), sizeof(int), 1, rep_data);

    fwrite(s->omega->row, sizeof(int), p->num->rv+1, rep_data);
    
    
    fwrite(s->omega->col, sizeof(int), p->num->rv+1, rep_data);
    
    for (idx = 0; idx < s->omega->cnt; idx++) {

        fwrite(s->omega->idx[idx], sizeof(int), p->num->cipher+1, rep_data);
        
    }

    fwrite(s->omega->weight, sizeof(int), s->omega->cnt, rep_data);
    

    fwrite(s->omega->filter, sizeof(int), s->omega->cnt, rep_data);

    
    fwrite(s->omega->RT+1, sizeof(double), sd_global->omegas.num_omega, rep_data);
    
    
    /* 2. Start storing lambda_type */
    fwrite(&(c->lambda->cnt), sizeof(int), 1, rep_data);

    fwrite(c->lambda->row, sizeof(int), p->num->rv_rows+1, rep_data);
    
    for (idx = 0; idx < c->lambda->cnt; idx++) {

        fwrite(c->lambda->val[idx], sizeof(double), p->num->rv_rows+1, rep_data);
        
    }
    
    /* 3. Start storing sigma_type */
    fwrite(&(c->sigma->cnt), sizeof(int), 1, rep_data);

    fwrite(c->sigma->col, sizeof(int), p->num->nz_cols+1, rep_data);
    
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fwrite(&(c->sigma->val[idx].R), sizeof(double), 1, rep_data);

        fwrite(c->sigma->val[idx].T, sizeof(double), p->num->nz_cols+1, rep_data);
        
    }

    fwrite(c->sigma->ck, sizeof(int), c->sigma->cnt, rep_data);
    

    fwrite(c->sigma->lamb, sizeof(int), c->sigma->cnt, rep_data);
    
    
    /* 4. Start storing theta_type */
    fwrite(&(c->theta->cnt), sizeof(int), 1, rep_data);

    fwrite(c->theta->last, sizeof(int), c->theta->cnt, rep_data);
    

    fwrite(c->theta->k, sizeof(double), c->theta->cnt, rep_data);
    

    fwrite(c->theta->p, sizeof(double), c->theta->cnt, rep_data);
    
    
    /* 5. Start storing delta_type */

    fwrite(s->delta->col, sizeof(int), p->num->rv_cols+1, rep_data);
    
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                fwrite(&(s->delta->val[idx][obs].R), sizeof(double), 1, rep_data);

                fwrite(s->delta->val[idx][obs].T, sizeof(double), p->num->rv_cols+1, rep_data);
                
            }
        }
    }
    
    /* 6. Start storing cuts */
    fwrite(&(c->cuts->cnt), sizeof(int), 1, rep_data);
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        fwrite(&(c->cuts->val[idx]->omega_cnt), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->cut_obs), sizeof(int), 1, rep_data);

        fwrite(c->cuts->val[idx]->istar, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
        
        fwrite(&(c->cuts->val[idx]->slack_cnt), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->cell_num), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->row_num), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->alpha), sizeof(double), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->alpha_incumb), sizeof(double), 1, rep_data);

        fwrite(c->cuts->val[idx]->beta, sizeof(double), p->num->mast_cols+1, rep_data);
        
        fwrite(&(c->cuts->val[idx]->subfeaflag), sizeof(BOOL), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->is_incumbent), sizeof(BOOL), 1, rep_data);
        if (c->cuts->val[idx]->is_incumbent) {

            fwrite(c->cuts->val[idx]->subobj_omega, sizeof(double), c->cuts->val[idx]->omega_cnt, rep_data);
            

            fwrite(c->cuts->val[idx]->subobj_freq, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
            
        }
        
    }
    
    /* 7. Start storing run_time */
    fwrite(&(s->run_time->total_time), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->iteration_time), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_master_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_subprob_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->full_test_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->argmax_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->iteration_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_master_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_subprob_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->full_test_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->argmax_accum), sizeof(double), 1, rep_data);
    
    
    
    /* 8. Start storing the remaining of cell structure */
    fwrite(&(c->id_num), sizeof(int), 1, rep_data);
    fwrite(&(c->num_members), sizeof(int), 1, rep_data);

    fwrite(c->members, sizeof(int), c->num_members, rep_data);
    
    fwrite(&(c->quad_scalar), sizeof(double), 1, rep_data);
    fwrite(&(c->LP_cnt), sizeof(int), 1, rep_data);
    fwrite(&(c->LP_test), sizeof(int), 1, rep_data);
    fwrite(&(c->N), sizeof(int), 1, rep_data);
    fwrite(&(c->P), sizeof(double), 1, rep_data);
    fwrite(&(c->k), sizeof(int), 1, rep_data);
    fwrite(&(c->opt_mode), sizeof(BOOL), 1, rep_data);
    fwrite(&(c->incumb_infea), sizeof(BOOL), 1, rep_data);
    fwrite(&(c->fea_count), sizeof(int), 1, rep_data);
    
    /* 9. Start storing the remaining of soln structure */

    fwrite(s->Pi, sizeof(double), p->num->sub_rows+1, rep_data);
    
    fwrite(&(s->subobj_est), sizeof(double), 1, rep_data);

    fwrite(s->Master_pi, sizeof(double), p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt + 1, rep_data);
    

    fwrite(s->Master_dj, sizeof(double), p->num->mast_cols + 1, rep_data);
    

    fwrite(s->candid_x, sizeof(double), p->num->mast_cols + 1, rep_data);

    fwrite(&(s->candid_est), sizeof(double), 1, rep_data);

    fwrite(s->incumb_x, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fwrite(&(s->incumb_k), sizeof(int), 1, rep_data);

    fwrite(s->incumb_d, sizeof(double), p->num->mast_cols + 1, rep_data);
    

    fwrite(s->incumb_avg, sizeof(double), p->num->mast_cols + 1, rep_data);
 
    fwrite(&(s->alpha), sizeof(double), 1, rep_data);

    fwrite(s->beta, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fwrite(&(s->incumb_est), sizeof(double), 1, rep_data);
    fwrite(&(s->opt_value), sizeof(double), 1, rep_data);
    fwrite(&(s->norm_d_k_1), sizeof(double), 1, rep_data);
    fwrite(&(s->norm_d_k), sizeof(double), 1, rep_data);
    fwrite(&(s->incumb_stdev), sizeof(double), 1, rep_data);
    fwrite(&(s->incumb_cut), sizeof(int), 1, rep_data);
    fwrite(&(s->last_update), sizeof(int), 1, rep_data);
    fwrite(&(s->gamma), sizeof(double), 1, rep_data);
    fwrite(&(s->optimality_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->smpl_ever_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->smpl_test_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->incumbent_change), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->dual_statble_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->full_test_error), sizeof(double), 1, rep_data);
    fwrite(&(s->passed), sizeof(int), 1, rep_data);
    fwrite(&(s->sub_lb_checker), sizeof(double), 1, rep_data);
    fwrite(&(s->max_ratio), sizeof(double), 1, rep_data);
    fwrite(&(s->min_ratio), sizeof(double), 1, rep_data);

    fwrite(s->pi_ratio, sizeof(double), sd_global->config.SCAN_LEN, rep_data);
    
    
    /* 10. Start storing the last seed used in getting observation */
    fwrite(&(sd_global->config.RUN_SEED), sizeof(sd_long), 1, rep_data);
    fwrite(&(sd_global->pi_flag[0]), sizeof(BOOL), 1, rep_data);
    fwrite(&(sd_global->pi_flag[1]), sizeof(BOOL), 1, rep_data);
    fwrite(&(sd_global->pi_flag[2]), sizeof(BOOL), 1, rep_data);
    
    /* 11. Start stroing the whole master problem */
    /* get the objective coefficient */
    for (idx = 0; idx <= p->master->mac; idx++) {
        get_coef(c->master, -1, idx, &coef_value);
        fwrite(&coef_value, sizeof(double), 1, rep_data);
    }
    
    /* get the first stage origianl constraint coefficient */
    for (cnt = 0; cnt < p->master->mar; cnt++) {
        for (idx = 0; idx < p->master->mac; idx++) {
            get_coef(c->master, cnt, idx, &coef_value);
            fwrite(&coef_value, sizeof(double), 1, rep_data);
        }
    }
    
    /* get the first stage origianl constraint rhs */
    for (cnt = 0; cnt < p->master->mar; cnt++) {
        get_coef(c->master, cnt, -1, &coef_value);
        fwrite(&coef_value, sizeof(double), 1, rep_data);
    }
    
    /* get the first stage added cuts' coefficients */
    for (cnt = p->master->mar; cnt < p->master->mar + c->cuts->cnt + c->feasible_cuts_added->cnt; cnt++) {
        /* Do not foget to inlcude the eta column*/
        for (idx = 0; idx <= p->master->mac; idx++) {
            get_coef(c->master, cnt, idx, &coef_value);
            fwrite(&coef_value, sizeof(double), 1, rep_data);
        }
        
    }
    
    /* get the first stage added cuts' rhs */
    for (cnt = p->master->mar; cnt < p->master->mar + c->cuts->cnt + c->feasible_cuts_added->cnt; cnt++) {
        get_coef(c->master, cnt, -1, &coef_value);
        fwrite(&coef_value, sizeof(double), 1, rep_data);
    }
    
    /* get the first stage decisions' upper bounds */
    for (idx = 0; idx < p->master->mac; idx++) {
        get_ubound(c->master, &coef_value, idx, idx);
        fwrite(&coef_value, sizeof(double), 1, rep_data);
    }
    
    /* get the first stage decisions' lower bounds */
    for (idx = 0; idx < p->master->mac; idx++) {
        get_lbound(c->master, &coef_value, idx, idx);
        fwrite(&coef_value, sizeof(double), 1, rep_data);
    }

    
    fclose(rep_data);
    
    /* 12. Set store_flag to TRUE indicating that no more stroring during this replication */
    /* This will be reset to FALSE at the beginning of each replication */
    sd_global->store_flag = TRUE;
    
    return 0;
}

int restore_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int cnt, idx, obs;
    int omega_cnt, cut_obs;
    FILE *rep_data;
    double coef_value;
    
    char buffer1[128];
//    char buffer2[128];
    char rep_number[16];
    sprintf(rep_number, "%d", p->current_batch_id);
    // sprintf(rep_number, "%d", 14);
    
    /* modified by Yifan 2014.02.26 Numbering resume data for different replications */
//    strcpy(buffer2, "./sdresume/pgp2/nominal/resume");
//    strcat(buffer2, rep_number);
//    strcat(buffer2, ".lp");
//    read_problem_simple(c->master, buffer2, "lp");
    
    /* modified by Yifan 2014.02.26 */
    strcpy(buffer1, "./sdresume/pgp2/nominal/resume_data");
    strcat(buffer1, rep_number);
    strcat(buffer1, ".txt");
    
    
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen(buffer1, "r");
    
    /* 1. Start restoring omega_type */
    
    if (fread(&(s->omega->cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->omega->cnt");
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        s->omega->idx[idx] = arr_alloc(p->num->cipher+1, int);
    }
    if (fread(&(s->omega->next), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->omega->next");
    }
    if ( fread(&(s->omega->most), sizeof(int), 1, rep_data)!= 1) {
        printf("Failed to read s->omega->most");
    }
    if ( fread(&(s->omega->last), sizeof(int), 1, rep_data)!= 1) {
        printf("Failed to read s->omega->last");
    }
    if ( fread(&(s->omega->k), sizeof(int), 1, rep_data)!= 1) {
        printf("Failed to read s->omega->k");
    }
    
    if (fread(s->omega->row, sizeof(int), p->num->rv+1, rep_data) != p->num->rv+1) {
        printf("Failed to read s->omega->row");
    }
    if (fread(s->omega->col, sizeof(int), p->num->rv+1, rep_data) != p->num->rv+1) {
        printf("Failed to read s->omega->col");
    }
    
    for (idx = 0; idx < s->omega->cnt; idx++) {
        
        if (fread(s->omega->idx[idx], sizeof(int), p->num->cipher+1, rep_data) != p->num->cipher+1) {
            printf("Failed to read s->omega->idx");
        }
        
    }
    
    if (fread(s->omega->weight, sizeof(int), s->omega->cnt, rep_data) != s->omega->cnt) {
        printf("Failed to read s->omega->weight");
    }
    
    
    if (fread(s->omega->filter, sizeof(int), s->omega->cnt, rep_data) != s->omega->cnt) {
        printf("Failed to read s->omega->filter");
    }
    
    
    if (fread(s->omega->RT+1, sizeof(double), sd_global->omegas.num_omega, rep_data) != sd_global->omegas.num_omega) {
        printf("Failed to read s->omega->RT+1");
    }
    
    
    /* 2. Start restoring lambda_type */
    if (fread(&(c->lambda->cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->lambda->cnt");
    }
    
    /* Dynamically allocate memory space for lambda */
    for (cnt = 0; cnt < c->lambda->cnt; cnt++)
		if (!(c->lambda->val[cnt] = arr_alloc(p->num->rv_rows+1, double)))
			err_msg("Allocation", "restore_lambda", "lambda->val[cnt]");
    

    if (fread(c->lambda->row, sizeof(int), p->num->rv_rows+1, rep_data) != p->num->rv_rows+1) {
        printf("Failed to read c->lambda->row");
    }

    for (idx = 0; idx < c->lambda->cnt; idx++) {
        
        if (fread(c->lambda->val[idx], sizeof(double), p->num->rv_rows+1, rep_data) != p->num->rv_rows+1) {
            printf("Failed to read c->lambda->val");
        }
        
    }
    
    /* 3. Start restoring sigma_type */
    if (fread(&(c->sigma->cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->sigma->cnt");
    }
    
    /* Dynamically allocate memory space for sigma */
    for (cnt = 0; cnt < c->sigma->cnt && cnt < sd_global->config.MAX_ITER; cnt++)
		if (!(c->sigma->val[cnt].T = arr_alloc(p->num->nz_cols+1, double)))
			err_msg("Allocation", "resotre_sigma", "sigma->val[cnt]");
    
    
    if (fread(c->sigma->col, sizeof(int), p->num->nz_cols+1, rep_data) != p->num->nz_cols+1) {
        printf("Failed to read c->sigma->col");
    }

    for (idx = 0; idx < c->sigma->cnt; idx++) {
        if (fread(&(c->sigma->val[idx].R), sizeof(double), 1, rep_data) != 1) {
            printf("Failed to read c->sigma->val[idx].R");
        }
        
        if (fread(c->sigma->val[idx].T, sizeof(double), p->num->nz_cols+1, rep_data) != p->num->nz_cols+1) {
            printf("Failed to read c->sigma->val[idx].T");
        }
        
    }
    
    if (fread(c->sigma->ck, sizeof(int), c->sigma->cnt, rep_data) != c->sigma->cnt) {
        printf("Failed to read c->sigma->ck");
    }

    
    if (fread(c->sigma->lamb, sizeof(int), c->sigma->cnt, rep_data) != c->sigma->cnt) {
        printf("Failed to read c->sigma->lamb");
    }

    
    /* 4. Start restoring theta_type */
    if (fread(&(c->theta->cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->theta->cnt");
    }
    
    if (fread(c->theta->last, sizeof(int), c->theta->cnt, rep_data) != c->theta->cnt) {
        printf("Failed to read c->theta->last");
    }

    
    if (fread(c->theta->k, sizeof(double), c->theta->cnt, rep_data) != c->theta->cnt) {
        printf("Failed to read c->theta->k");
    }

    
    if (fread(c->theta->p, sizeof(double), c->theta->cnt, rep_data) != c->theta->cnt) {
        printf("Failed to read c->theta->p");
    }

    
    
    /* 5. Start restoring delta_type */
    
    /* Dynamically allocate memory space for delta */
    /* Double check this part since there might be an error here!!!! */
    for (cnt = 0; cnt < c->lambda->cnt; cnt++){
        if (!(s->delta->val[cnt] = arr_alloc(p->num->iter, pi_R_T_type)))
            err_msg("Allocation", "restore_delta", "delta->val");
        for (obs=0; obs<s->omega->most; obs++) {
            if (!(s->delta->val[cnt][obs].T = arr_alloc(p->num->rv_cols + 1, double)))
                err_msg("Allocation", "restore_delta", "delta->val[cnt]");
        }
    }
    
    
    
    if (fread(s->delta->col, sizeof(int), p->num->rv_cols+1, rep_data) != p->num->rv_cols+1) {
        printf("Failed to read s->delta->col");
    }

    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                if (fread(&(s->delta->val[idx][obs].R), sizeof(double), 1, rep_data) != 1) {
                    printf("Failed to read s->delta->val[idx][obs].R");
                }
                
                if (fread(s->delta->val[idx][obs].T, sizeof(double), p->num->rv_cols+1, rep_data) != p->num->rv_cols+1) {
                    printf("Failed to read s->delta->val[idx][obs].T");
                }
                
            }
        }
    }
    
    /* 6. Start restoring cuts */
    if (fread(&(c->cuts->cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->cuts->cnt");
    }
    
    /* Dynamically allocate memory space for cuts */
    
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        if (fread(&(omega_cnt), sizeof(int), 1, rep_data) != 1) {
            printf("Failed to read omega_cnt");
        }
        if (fread(&(cut_obs), sizeof(int), 1, rep_data) != 1) {
            printf("Failed to read cut_obs");
        }
        c->cuts->val[idx] = new_cut(p->num->mast_cols, omega_cnt, cut_obs);
        c->cuts->val[idx]->omega_cnt = omega_cnt;
        c->cuts->val[idx]->cut_obs = cut_obs;
        
        if (fread(c->cuts->val[idx]->istar, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data) != c->cuts->val[idx]->omega_cnt) {
            printf("Failed to read c->cuts->val[idx]->istar");
        }
        
        if (fread(&(c->cuts->val[idx]->slack_cnt), sizeof(int), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->slack_cnt");
        }
        if (fread(&(c->cuts->val[idx]->cell_num), sizeof(int), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->cell_num");
        }
        if (fread(&(c->cuts->val[idx]->row_num), sizeof(int), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->row_num");
        }
        if (fread(&(c->cuts->val[idx]->alpha), sizeof(double), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->alpha");
        }
        if (fread(&(c->cuts->val[idx]->alpha_incumb), sizeof(double), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->alpha_incumb");
        }
        
        if (fread(c->cuts->val[idx]->beta, sizeof(double), p->num->mast_cols+1, rep_data) != p->num->mast_cols+1) {
            printf("Failed to read c->cuts->val[idx]->beta");
        }
        
        if (fread(&(c->cuts->val[idx]->subfeaflag), sizeof(BOOL), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->subfeaflag");
        }
        if (fread(&(c->cuts->val[idx]->is_incumbent), sizeof(BOOL), 1, rep_data) != 1) {
            printf("Failed to read c->cuts->val[idx]->is_incumbent");
        }
        if (c->cuts->val[idx]->is_incumbent) {
            if (!(c->cuts->val[idx]->subobj_omega = arr_alloc( c->cuts->val[idx]->omega_cnt+1, double)))
                err_msg("Allocation", "SD_cut", "subobj_omega");
            if (!(c->cuts->val[idx]->subobj_freq = arr_alloc( c->cuts->val[idx]->omega_cnt+1, int)))
                err_msg("Allocation", "SD_cut", "subobj_freq");
            
            if (fread(c->cuts->val[idx]->subobj_omega, sizeof(double), c->cuts->val[idx]->omega_cnt, rep_data) != c->cuts->val[idx]->omega_cnt) {
                printf("Failed to read c->cuts->val[idx]->subobj_omega");
            }
            
            
            if (fread(c->cuts->val[idx]->subobj_freq, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data) != c->cuts->val[idx]->omega_cnt) {
                printf("Failed to read c->cuts->val[idx]->subobj_freq");
            }
            
        }
        
    }
    
    /* 7. Start restoring run_time */
    if (fread(&(s->run_time->total_time), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->total_time");
    }
    if (fread(&(s->run_time->iteration_time), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->iteration_time");
    }
    if (fread(&(s->run_time->soln_master_iter), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->soln_master_iter");
    }
    if (fread(&(s->run_time->soln_subprob_iter), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->soln_subprob_iter");
    }
    if (fread(&(s->run_time->full_test_iter), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->full_test_iter");
    }
    if (fread(&(s->run_time->argmax_iter), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->argmax_iter");
    }
    if (fread(&(s->run_time->iteration_accum), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->iteration_accum");
    }
    if (fread(&(s->run_time->soln_master_accum), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->soln_master_accum");
    }
    if (fread(&(s->run_time->soln_subprob_accum), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->soln_subprob_accum");
    }
    if (fread(&(s->run_time->full_test_accum), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->full_test_accum");
    }
    if (fread(&(s->run_time->argmax_accum), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->run_time->argmax_accum");
    }
    
    
    
    /* 8. Start restoring the remaining of cell structure */
    if (fread(&(c->id_num), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->id_num");
    }
    if (fread(&(c->num_members), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->num_members");
    }
    
    if (fread(c->members, sizeof(int), c->num_members, rep_data) != c->num_members) {
        printf("Failed to read c->members");
    }
    
    if (fread(&(c->quad_scalar), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read c->quad_scalar");
    }
    if (fread(&(c->LP_cnt), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->LP_cnt");
    }
    if (fread(&(c->LP_test), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->LP_test");
    }
    if (fread(&(c->N), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->N");
    }
    if (fread(&(c->P), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read c->P");
    }
    if (fread(&(c->k), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->k");
    }
    if (fread(&(c->opt_mode), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read c->opt_mode");
    }
    if (fread(&(c->incumb_infea), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read c->incumb_infea");
    }
    if (fread(&(c->fea_count), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read c->fea_count");
    }
    
    /* 9. Start restoring the remaining of soln structure */
    
    if (fread(s->Pi, sizeof(double), p->num->sub_rows+1, rep_data) != p->num->sub_rows+1) {
        printf("Failed to read s->Pi");
    }
    
    if (fread(&(s->subobj_est), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->subobj_est");
    }
    /* Update the dual size before restoring added by Yifan 2014.01.28*/
    update_dual_size(c, s, p);
    
    if (fread(s->Master_pi, sizeof(double), p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt + 1, rep_data) != p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt + 1) {
        printf("Failed to read s->Master_pi");
    }
    

    
    if (fread(s->Master_dj, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->Master_dj");
    }
    
    
    if (fread(s->candid_x, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->candid_x");
    }
    
    if (fread(&(s->candid_est), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->candid_est");
    }
    
    if (fread(s->incumb_x, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->incumb_x");
    }
    
    if (fread(&(s->incumb_k), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->incumb_k");
    }
    
    if (fread(s->incumb_d, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->incumb_d");
    }
    
    
    if (fread(s->incumb_avg, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->incumb_avg");
    }
    
    if (fread(&(s->alpha), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->alpha");
    }
    
    if (fread(s->beta, sizeof(double), p->num->mast_cols + 1, rep_data) != p->num->mast_cols + 1) {
        printf("Failed to read s->beta");
    }
    
    if (fread(&(s->incumb_est), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->incumb_est");
    }
    if (fread(&(s->opt_value), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->opt_value");
    }
    if (fread(&(s->norm_d_k_1), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->norm_d_k_1");
    }
    if (fread(&(s->norm_d_k), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->norm_d_k");
    }
    if (fread(&(s->incumb_stdev), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->incumb_stdev");
    }
    if (fread(&(s->incumb_cut), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->incumb_cut");
    }
    if (fread(&(s->last_update), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->last_update");
    }
    if (fread(&(s->gamma), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->gamma");
    }
    if (fread(&(s->optimality_flag), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read s->optimality_flag");
    }
    if (fread(&(s->smpl_ever_flag), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read s->smpl_ever_flag");
    }
    if (fread(&(s->smpl_test_flag), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read s->smpl_test_flag");
    }
    if (fread(&(s->incumbent_change), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read s->incumbent_change");
    }
    if (fread((s->dual_statble_flag), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read s->dual_statble_flag");
    }
    /* make sure dual_stable_flag is turned to NO */
    *s->dual_statble_flag = 0;
    if (fread(&(s->full_test_error), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->full_test_error");
    }
    if (fread(&(s->passed), sizeof(int), 1, rep_data) != 1) {
        printf("Failed to read s->passed");
    }
    if (fread(&(s->sub_lb_checker), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->sub_lb_checker");
    }
    if (fread(&(s->max_ratio), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->max_ratio");
    }
    if (fread(&(s->min_ratio), sizeof(double), 1, rep_data) != 1) {
        printf("Failed to read s->min_ratio");
    }
    
    if (fread(s->pi_ratio, sizeof(double), sd_global->config.SCAN_LEN, rep_data) != sd_global->config.SCAN_LEN) {
        printf("Failed to read s->pi_ratio");
    }
    
    
    /* 10. Start restoring the last seed used in getting observation */
    if (fread(&(sd_global->config.RUN_SEED), sizeof(sd_long), 1, rep_data) != 1) {
        printf("Failed to read sd_global->config.RUN_SEED");
    }
    if (fread(&(sd_global->pi_flag[0]), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read sd_global->pi_flag[0]");
    }
    if (fread(&(sd_global->pi_flag[1]), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read sd_global->pi_flag[1]");
    }
    if (fread(&(sd_global->pi_flag[2]), sizeof(BOOL), 1, rep_data) != 1) {
        printf("Failed to read sd_global->pi_flag[2]");
    }
    
    /* 11. Start stroing the whole master problem */
    if (0) {
        /* get the objective coefficient */
        for (idx = 0; idx <= p->master->mac; idx++) {
            if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                printf("Failed to read obj coef_value");
            }
            change_single_coef(c->master, -1, idx, coef_value);
        }
        
        /* get the first stage origianl constraint coefficient */
        for (cnt = 0; cnt < p->master->mar; cnt++) {
            for (idx = 0; idx < p->master->mac; idx++) {
                if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                    printf("Failed to read constraint coef_value");
                }
                change_single_coef(c->master, cnt, idx, coef_value);
            }
        }
        
        /* get the first stage origianl constraint rhs */
        for (cnt = 0; cnt < p->master->mar; cnt++) {
            if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                printf("Failed to read rhs");
            }
            change_single_coef(c->master, cnt, -1, coef_value);
        }
        
        /* get the first stage added cuts' coefficients */
        for (cnt = p->master->mar; cnt < p->master->mar + c->cuts->cnt + c->feasible_cuts_added->cnt; cnt++) {
            /* Do not foget to inlcude the eta column*/
            for (idx = 0; idx <= p->master->mac; idx++) {
                if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                    printf("Failed to read eta coef_value");
                }
                change_single_coef(c->master, cnt, idx, coef_value);
            }
            
        }
        
        /* get the first stage added cuts' rhs */
        for (cnt = p->master->mar; cnt < p->master->mar + c->cuts->cnt + c->feasible_cuts_added->cnt; cnt++) {
            if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                printf("Failed to read cut rhs");
            }
            change_single_coef(c->master, cnt, -1, coef_value);
        }
        
        /* get the first stage decisions' upper bounds */
        for (idx = 0; idx < p->master->mac; idx++) {
            if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                printf("Failed to read x's upper bound");
            }
            change_bound(c->master, 1, &idx, "U", &coef_value);
        }
        
        /* get the first stage decisions' lower bounds */
        for (idx = 0; idx < p->master->mac; idx++) {
            if (fread(&coef_value, sizeof(double), 1, rep_data) != 1) {
                printf("Failed to read x's lower bound");
            }
            change_bound(c->master, 1, &idx, "L", &coef_value);
        }
    }

    
    fclose(rep_data);
    
    return 0;
}
