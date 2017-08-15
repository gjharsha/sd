//
//  resume.c
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
#include "resume.h"

int store_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int cnt, idx, obs;
    FILE *rep_data;
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen("rep_data.txt", "w");
    /* 1. Start storing omega_type */
    fprintf(rep_data, "1. Start storing omega_type\n");
    fprintf(rep_data, "%d\n",s->omega->cnt);
    fprintf(rep_data, "%d\n",s->omega->next);
    fprintf(rep_data, "%d\n",s->omega->most);
    fprintf(rep_data, "%d\n",s->omega->last);
    fprintf(rep_data, "%d\n",s->omega->k);
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        fprintf(rep_data, "%d\n", s->omega->row[cnt]);
    }
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        fprintf(rep_data, "%d\n", s->omega->col[cnt]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->cipher; cnt++) {
            fprintf(rep_data, "%d\n", s->omega->idx[idx][cnt]);
        }
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        fprintf(rep_data, "%d\n", s->omega->weight[idx]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        fprintf(rep_data, "%d\n", s->omega->filter[idx]);
    }
    for (cnt = 1; cnt <= sd_global->omegas.num_omega; cnt++) {
        fprintf(rep_data, "%.17g\n", s->omega->RT[cnt]);
    }
    
    /* 2. Start storing lambda_type */
    fprintf(rep_data, "2. Start storing lambda_type\n");
    fprintf(rep_data, "%d\n",c->lambda->cnt);
    for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
        fprintf(rep_data, "%d\n", c->lambda->row[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
            fprintf(rep_data, "%.17g\n", c->lambda->val[idx][cnt]);
        }
    }
    
    /* 3. Start storing sigma_type */
    fprintf(rep_data, "3. Start storing sigma_type\n");
    fprintf(rep_data, "%d\n",c->sigma->cnt);
    for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
        fprintf(rep_data, "%d\n", c->sigma->col[cnt]);
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fprintf(rep_data, "%.17g\n", c->sigma->val[idx].R);
        for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
            fprintf(rep_data, "%.17g\n", c->sigma->val[idx].T[cnt]);
        }
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fprintf(rep_data, "%d\n", c->sigma->ck[idx]);
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fprintf(rep_data, "%d\n", c->sigma->lamb[idx]);
    }
    
    /* 4. Start storing theta_type */
    fprintf(rep_data, "4. Start storing theta_type\n");
    fprintf(rep_data, "cnt: %d\n",c->theta->cnt);
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fprintf(rep_data, "%d\n", c->theta->last[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fprintf(rep_data, "%.17g\n", c->theta->k[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fprintf(rep_data, "%.17g\n", c->theta->p[idx]);
    }
    
    /* 5. Start storing delta_type */
    fprintf(rep_data, "5. Start storing delta_type\n");
    for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
        fprintf(rep_data, "%d\n", s->delta->col[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                fprintf(rep_data, "%.17g\n", s->delta->val[idx][obs].R);
                for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
                    fprintf(rep_data, "%.17g\n", s->delta->val[idx][obs].T[cnt]);
                }
            }
        }
    }

    /* 6. Start storing cuts */
    fprintf(rep_data, "6. Start storing cuts\n");
    fprintf(rep_data, "cnt: %d\n", c->cuts->cnt);
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        fprintf(rep_data, "%d\n", c->cuts->val[idx]->omega_cnt);
        fprintf(rep_data, "%d\n", c->cuts->val[idx]->cut_obs);
        for (cnt = 0; cnt < c->cuts->val[idx]->omega_cnt; cnt++) {
            fprintf(rep_data, "%d\n", c->cuts->val[idx]->istar[cnt]);
        }
        fprintf(rep_data, "%d\n", c->cuts->val[idx]->slack_cnt);
        fprintf(rep_data, "%d\n", c->cuts->val[idx]->cell_num);
        fprintf(rep_data, "%d\n", c->cuts->val[idx]->row_num);
        fprintf(rep_data, "%.17g\n", c->cuts->val[idx]->alpha);
        fprintf(rep_data, "%.17g\n", c->cuts->val[idx]->alpha_incumb);
        for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
            fprintf(rep_data, "%.17g\n", c->cuts->val[idx]->beta[cnt]);
        }
        fprintf(rep_data, "%u\n", c->cuts->val[idx]->subfeaflag);
        fprintf(rep_data, "%u\n", c->cuts->val[idx]->is_incumbent);
        if (c->cuts->val[idx]->is_incumbent) {
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                fprintf(rep_data, "%.17g\n", c->cuts->val[idx]->subobj_omega[obs]);
            }
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                fprintf(rep_data, "%d\n", c->cuts->val[idx]->subobj_freq[obs]);
            }
        }
        
    }
    
    /* 7. Start storing run_time */
    fprintf(rep_data, "7. Start storing run_time\n");
    fprintf(rep_data, "%.17g\n", s->run_time->total_time);
    fprintf(rep_data, "%.17g\n", s->run_time->iteration_time);
    fprintf(rep_data, "%.17g\n", s->run_time->soln_master_iter);
    fprintf(rep_data, "%.17g\n", s->run_time->soln_subprob_iter);
    fprintf(rep_data, "%.17g\n", s->run_time->full_test_iter);
    fprintf(rep_data, "%.17g\n", s->run_time->argmax_iter);
    fprintf(rep_data, "%.17g\n", s->run_time->iteration_accum);
    fprintf(rep_data, "%.17g\n", s->run_time->soln_master_accum);
    fprintf(rep_data, "%.17g\n", s->run_time->soln_subprob_accum);
    fprintf(rep_data, "%.17g\n", s->run_time->full_test_accum);
    fprintf(rep_data, "%.17g\n", s->run_time->argmax_accum);



    /* 8. Start storing the remaining of cell structure */
    fprintf(rep_data, "8. Start storing the remaining of cell structure\n");
    fprintf(rep_data, "%d\n", c->id_num);
    fprintf(rep_data, "%d\n", c->num_members);
    for (idx = 0; idx < c->num_members; idx++) {
        fprintf(rep_data, "%d\n", c->members[idx]);
    }
    fprintf(rep_data, "%.17g\n", c->quad_scalar);
    fprintf(rep_data, "%d\n", c->LP_cnt);
    fprintf(rep_data, "%d\n", c->LP_test);
    fprintf(rep_data, "%d\n", c->N);
    fprintf(rep_data, "%.17g\n", c->P);
    fprintf(rep_data, "%d\n", c->k);
    fprintf(rep_data, "%u\n", c->opt_mode);
    fprintf(rep_data, "%u\n", c->incumb_infea);
    fprintf(rep_data, "%d\n", c->fea_count);
    
    /* 9. Start storing the remaining of soln structure */
    fprintf(rep_data, "9. Start storing the remaining of soln structure\n");
    fprintf(rep_data, "9-1. Start storing s->Pi\n");
    for (idx = 0; idx <= p->num->sub_rows; idx++) {
        fprintf(rep_data, "%.17g\n", s->Pi[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->subobj_est\n");
    fprintf(rep_data, "%.17g\n", s->subobj_est);
    fprintf(rep_data, "9-1. Start storing s->Master_pi\n");
    for (idx = 0; idx <= p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt; idx++) {
        fprintf(rep_data, "%.17g\n", s->Master_pi[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->Master_dj\n");
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fprintf(rep_data, "%.17g\n", s->Master_dj[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->candid_x\n");
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fprintf(rep_data, "%.17g\n", s->candid_x[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->candid_est\n");
    fprintf(rep_data, "%.17g\n", s->candid_est);
    fprintf(rep_data, "9-1. Start storing s->incumb_x\n");
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fprintf(rep_data, "%.17g\n", s->incumb_x[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->incumb_k\n");
    fprintf(rep_data, "%d\n", s->incumb_k);
    fprintf(rep_data, "9-1. Start storing s->incumb_d\n");
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fprintf(rep_data, "%.17g\n", s->incumb_d[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->incumb_avg\n");
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fprintf(rep_data, "%.17g\n", s->incumb_avg[idx]);
    }
    fprintf(rep_data, "9-1. Start storing s->alpha\n");
    fprintf(rep_data, "%.17g\n", s->alpha);
    fprintf(rep_data, "9-1. Start storing s->beta\n");
    for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
        fprintf(rep_data, "%.17g\n", s->beta[cnt]);
    }
    fprintf(rep_data, "9-1. Start storing single value\n");
    fprintf(rep_data, "%.17g\n", s->incumb_est);
    fprintf(rep_data, "%.17g\n", s->opt_value);
    fprintf(rep_data, "%.17g\n", s->norm_d_k_1);
    fprintf(rep_data, "%.17g\n", s->norm_d_k);
    fprintf(rep_data, "%.17g\n", s->incumb_stdev);
    fprintf(rep_data, "%d\n", s->incumb_cut);
    fprintf(rep_data, "%d\n", s->last_update);
    fprintf(rep_data, "%.17g\n", s->gamma);
    fprintf(rep_data, "%u\n", s->optimality_flag);
    fprintf(rep_data, "%u\n", s->smpl_ever_flag);
    fprintf(rep_data, "%u\n", s->smpl_test_flag);
    fprintf(rep_data, "%u\n", s->incumbent_change);
    fprintf(rep_data, "%u\n", *(s->dual_statble_flag));
    fprintf(rep_data, "%.17g\n", s->full_test_error);
    fprintf(rep_data, "%d\n", s->passed);
    fprintf(rep_data, "%.17g\n", s->sub_lb_checker);
    fprintf(rep_data, "%.17g\n", s->max_ratio);
    fprintf(rep_data, "%.17g\n", s->min_ratio);
    fprintf(rep_data, "9-1. Start storing s->pi_ratio\n");
    for (idx = 0; idx < sd_global->config.SCAN_LEN; idx++) {
        fprintf(rep_data, "%.17g\n", s->pi_ratio[idx]);
    }

    /* 10. Start storing the last seed used in getting observation */
    fprintf(rep_data, "10. Start storing the last seed used in getting observation\n");
    fprintf(rep_data, "%lld\n", sd_global->config.RUN_SEED);
    
    fclose(rep_data);
    
    return 0;
}

int restore_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int cnt, idx, obs;
    int omega_cnt, cut_obs;
    FILE *rep_data;
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen("rep_data.txt", "r");
    
    /* 1. Start restoring omega_type */
    fscanf(rep_data, "%d\n",&s->omega->cnt);
    for (idx = 0; idx < s->omega->cnt; idx++) {
        s->omega->idx[idx] = arr_alloc(p->num->cipher+1, int);
    }
    
    fscanf(rep_data, "%d\n",&s->omega->next);
    fscanf(rep_data, "%d\n",&s->omega->most);
    fscanf(rep_data, "%d\n",&s->omega->last);
    fscanf(rep_data, "%d\n",&s->omega->k);
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        fscanf(rep_data, "%d\n", &s->omega->row[cnt]);
    }
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        fscanf(rep_data, "%d\n", &s->omega->col[cnt]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->cipher; cnt++) {
            fscanf(rep_data, "%d\n", &s->omega->idx[idx][cnt]);
        }
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        fscanf(rep_data, "%d\n", &s->omega->weight[idx]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        fscanf(rep_data, "%d\n", &s->omega->filter[idx]);
    }
    for (cnt = 1; cnt <= sd_global->omegas.num_omega; cnt++) {
        fscanf(rep_data, "%lg\n", &s->omega->RT[cnt]);
    }
    
    /* 2. Start restoring lambda_type */
    fscanf(rep_data, "%d\n",&c->lambda->cnt);
    
    /* Dynamically allocate memory space for lambda */
    for (cnt = 0; cnt < c->lambda->cnt; cnt++)
		if (!(c->lambda->val[cnt] = arr_alloc(p->num->rv_rows+1, double)))
			err_msg("Allocation", "restore_lambda", "lambda->val[cnt]");
    
    for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
        fscanf(rep_data, "%d\n", &c->lambda->row[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
            fscanf(rep_data, "%lg\n", &c->lambda->val[idx][cnt]);
        }
    }
    
    /* 3. Start restoring sigma_type */
    fscanf(rep_data, "%d\n",&c->sigma->cnt);
    
    /* Dynamically allocate memory space for sigma */
    for (cnt = 0; cnt < c->sigma->cnt && cnt < sd_global->config.MAX_ITER; cnt++)
		if (!(c->sigma->val[cnt].T = arr_alloc(p->num->nz_cols+1, double)))
			err_msg("Allocation", "resotre_sigma", "sigma->val[cnt]");
    
    for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
        fscanf(rep_data, "%d\n", &c->sigma->col[cnt]);
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fscanf(rep_data, "%lg\n", &c->sigma->val[idx].R);
        for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
            fscanf(rep_data, "%lg\n", &c->sigma->val[idx].T[cnt]);
        }
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fscanf(rep_data, "%d\n", &c->sigma->ck[idx]);
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fscanf(rep_data, "%d\n", &c->sigma->lamb[idx]);
    }
    
    /* 4. Start restoring theta_type */
    fscanf(rep_data, "cnt: %d\n",&c->theta->cnt);
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fscanf(rep_data, "%d\n", &c->theta->last[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fscanf(rep_data, "%lg\n", &c->theta->k[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        fscanf(rep_data, "%lg\n", &c->theta->p[idx]);
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
    
    
    for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
        fscanf(rep_data, "%d\n", &s->delta->col[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                fscanf(rep_data, "%lg\n", &s->delta->val[idx][obs].R);
                for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
                    fscanf(rep_data, "%lg\n", &s->delta->val[idx][obs].T[cnt]);
                }
            }
        }
    }
    
    /* 6. Start restoring cuts */
    fscanf(rep_data, "cnt: %d\n", &c->cuts->cnt);
    
    /* Dynamically allocate memory space for cuts */
		
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        fscanf(rep_data, "%d\n", &omega_cnt);
        fscanf(rep_data, "%d\n", &cut_obs);
        c->cuts->val[idx] = new_cut(p->num->mast_cols, omega_cnt, cut_obs);
        c->cuts->val[idx]->omega_cnt = omega_cnt;
        c->cuts->val[idx]->cut_obs = cut_obs;
        for (cnt = 0; cnt < c->cuts->val[idx]->omega_cnt; cnt++) {
            fscanf(rep_data, "%d\n", &c->cuts->val[idx]->istar[cnt]);
        }
        fscanf(rep_data, "%d\n", &c->cuts->val[idx]->slack_cnt);
        fscanf(rep_data, "%d\n", &c->cuts->val[idx]->cell_num);
        fscanf(rep_data, "%d\n", &c->cuts->val[idx]->row_num);
        fscanf(rep_data, "%lg\n", &c->cuts->val[idx]->alpha);
        fscanf(rep_data, "%lg\n", &c->cuts->val[idx]->alpha_incumb);
        for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
            fscanf(rep_data, "%lg\n", &c->cuts->val[idx]->beta[cnt]);
        }
        fscanf(rep_data, "%u\n", &c->cuts->val[idx]->subfeaflag);
        fscanf(rep_data, "%u\n", &c->cuts->val[idx]->is_incumbent);
        if (c->cuts->val[idx]->is_incumbent) {
            if (!(c->cuts->val[idx]->subobj_omega = arr_alloc( c->cuts->val[idx]->omega_cnt+1, double)))
                err_msg("Allocation", "SD_cut", "subobj_omega");
            if (!(c->cuts->val[idx]->subobj_freq = arr_alloc( c->cuts->val[idx]->omega_cnt+1, int)))
                err_msg("Allocation", "SD_cut", "subobj_freq");
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                fscanf(rep_data, "%lg\n", &c->cuts->val[idx]->subobj_omega[obs]);
            }
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                fscanf(rep_data, "%d\n", &c->cuts->val[idx]->subobj_freq[obs]);
            }
        }
        
    }
    
    /* 6. Start restoring run_time */
    fscanf(rep_data, "%lg\n", &s->run_time->total_time);
    fscanf(rep_data, "%lg\n", &s->run_time->iteration_time);
    fscanf(rep_data, "%lg\n", &s->run_time->soln_master_iter);
    fscanf(rep_data, "%lg\n", &s->run_time->soln_subprob_iter);
    fscanf(rep_data, "%lg\n", &s->run_time->full_test_iter);
    fscanf(rep_data, "%lg\n", &s->run_time->argmax_iter);
    fscanf(rep_data, "%lg\n", &s->run_time->iteration_accum);
    fscanf(rep_data, "%lg\n", &s->run_time->soln_master_accum);
    fscanf(rep_data, "%lg\n", &s->run_time->soln_subprob_accum);
    fscanf(rep_data, "%lg\n", &s->run_time->full_test_accum);
    fscanf(rep_data, "%lg\n", &s->run_time->argmax_accum);
    
    
    
    /* 8. Start restoring the remaining of cell structure */
    fscanf(rep_data, "%d\n", &c->id_num);
    fscanf(rep_data, "%d\n", &c->num_members);
    for (idx = 0; idx < c->num_members; idx++) {
        fscanf(rep_data, "%d\n", &c->members[idx]);
    }
    fscanf(rep_data, "%lg\n", &c->quad_scalar);
    fscanf(rep_data, "%d\n", &c->LP_cnt);
    fscanf(rep_data, "%d\n", &c->LP_test);
    fscanf(rep_data, "%d\n", &c->N);
    fscanf(rep_data, "%lg\n", &c->P);
    fscanf(rep_data, "%d\n", &c->k);
    fscanf(rep_data, "%u\n", &c->opt_mode);
    fscanf(rep_data, "%u\n", &c->incumb_infea);
    fscanf(rep_data, "%d\n", &c->fea_count);
    
    /* 9. Start restoring the remaining of soln structure */
    for (idx = 0; idx <= p->num->sub_rows; idx++) {
        fscanf(rep_data, "%lg\n", &s->Pi[idx]);
    }
    fscanf(rep_data, "%lg\n", &s->subobj_est);
    /* Update the dual size before restoring added by Yifan 2014.01.28*/
    update_dual_size(c, s, p);
    for (idx = 0; idx <= p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt; idx++) {
        fscanf(rep_data, "%lg\n", &s->Master_pi[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fscanf(rep_data, "%lg\n", &s->Master_dj[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fscanf(rep_data, "%lg\n", &s->candid_x[idx]);
    }
    fscanf(rep_data, "%lg\n", &s->candid_est);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fscanf(rep_data, "%lg\n", &s->incumb_x[idx]);
    }
    fscanf(rep_data, "%d\n", &s->incumb_k);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fscanf(rep_data, "%lg\n", &s->incumb_d[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        fscanf(rep_data, "%lg\n", &s->incumb_avg[idx]);
    }
    fscanf(rep_data, "%lg\n", &s->alpha);
    for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
        fscanf(rep_data, "%lg\n", &s->beta[cnt]);
    }
    fscanf(rep_data, "%lg\n", &s->incumb_est);
    fscanf(rep_data, "%lg\n", &s->opt_value);
    fscanf(rep_data, "%lg\n", &s->norm_d_k_1);
    fscanf(rep_data, "%lg\n", &s->norm_d_k);
    fscanf(rep_data, "%lg\n", &s->incumb_stdev);
    fscanf(rep_data, "%d\n", &s->incumb_cut);
    fscanf(rep_data, "%d\n", &s->last_update);
    fscanf(rep_data, "%lg\n", &s->gamma);
    fscanf(rep_data, "%u\n", &s->optimality_flag);
    fscanf(rep_data, "%u\n", &s->smpl_ever_flag);
    fscanf(rep_data, "%u\n", &s->smpl_test_flag);
    fscanf(rep_data, "%u\n", &s->incumbent_change);
    fscanf(rep_data, "%u\n", (s->dual_statble_flag));
    /* make sure dual_stable_flag is turned to NO */
    *s->dual_statble_flag = 0;
    fscanf(rep_data, "%lg\n", &s->full_test_error);
    fscanf(rep_data, "%d\n", &s->passed);
    fscanf(rep_data, "%lg\n", &s->sub_lb_checker);
    fscanf(rep_data, "%lg\n", &s->max_ratio);
    fscanf(rep_data, "%lg\n", &s->min_ratio);
    for (idx = 0; idx < sd_global->config.SCAN_LEN; idx++) {
        fscanf(rep_data, "%lg\n", &s->pi_ratio[idx]);
    }
    
    /* 10. Start restoring the last seed used in getting observation */
    fscanf(rep_data, "%lld\n", &sd_global->config.RUN_SEED);
    
    fclose(rep_data);
    
    return 0;
}
