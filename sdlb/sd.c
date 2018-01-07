/***********************************************************************\
 ** sd.c
 ** 
 ** This program solves stochastic linear programs using stochastic
 ** decomposition.  The possible usages are:
 **
 **   sd	-- Prompts the user for one problem and solves it.
 **
 **   sd n -- Solves _n_ problems named "prob0000", "prob0001", etc.
 **
 **   sd n m -- Solves _n_ problems (as above) starting with "prob000m".
 **
 **   sd fname obj s1 s2 -- Solves the problem specified by _fname_ having
 **                         an objective sense of _obj_ (1 for minimization),
 **                         using random number seeds _s1_ and _s2_.
 **
 **
 ** This file contains the following functions:
 **   main()
 **   err_msg()
 **   parse_cmd_line()
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include <time.h>
#include "solver.h"
#include "utility.h"
#include "input.h"
#include "log.h"
#include "cuts.h"
#include "sdglobal.h"
#include "supomega.h"
#include "unistd.h"

#ifdef SD_win
#include <windows.h>
#endif

void sd_create_resume_folder(sdglobal_type* sd_global, char *buffer3, char *buffer2, char *fname);
BOOL sd_check_resume_folder(sdglobal_type* sd_global, char *fname);
BOOL sd_check_resume_data(char *buffer2, char *tolerance);
void sd_create_output_folder(sdglobal_type* sd_global, char *buffer1, char *buffer2, char *fname);
sdglobal_type *setupSDglobal();
void freeSDglobal(sdglobal_type *sd_global);

#ifdef SD_win
void sd_mv_output_files(char *fname);
#else
void sd_mv_resume_files(char *buffer1, char *buffer2, char *fname);
void sd_mv_output_files(char *buffer1, char *buffer2, char *fname);
#endif

/* Begin SD */
int main(int argc, char *argv[]) {
	sdglobal_type *SDglobal = NULL;
	one_problem *probptr = NULL;
	double *x_k = NULL, *original_x_k = NULL, objective;
	int num_rv, num_cipher, num_probs, n, start, objsen, cnt = 0, idx = 0;
	sd_small row, col;
	identity *ident = NULL;
	char fname[NAME_SIZE], buffer1[128], buffer2[128],buffer3[128];
	BOOL read_seeds, read_iters;
	sd_long seed1[BATCH_SIZE], seed2[BATCH_SIZE];

	/* Initialize the CPLEX environment, and to start the algorithm */
	SDglobal = setupSDglobal();
	open_Solver();
	log_start(SDglobal);

	/* parse the command line to read the input */
	parse_cmd_line(SDglobal, argc, argv, fname, &objsen, &num_probs, &start, &read_seeds, &read_iters);

	/* Load the solution settings */
	if (!load_config(SDglobal, read_seeds, read_iters))
		return 1;

#ifdef SD_unix
	if (RESUME_FLAG) {
		SDglobal->resume_flag = sd_check_resume_folder(SDglobal, fname);
		sd_create_resume_folder(SDglobal, buffer3,buffer2,fname);
	}
	sd_create_output_folder(SDglobal, buffer1,buffer2,fname);
#endif

	if (SDglobal->config.MULTIPLE_REP == 1) {
		if (SDglobal->config.AUTO_SEED == 0) {
			seed1[0] = SDglobal->config.RUN_SEED1;
			seed1[1] = SDglobal->config.RUN_SEED2;
			seed1[2] = SDglobal->config.RUN_SEED3;
			seed1[3] = SDglobal->config.RUN_SEED4;
			seed1[4] = SDglobal->config.RUN_SEED5;
			seed1[5] = SDglobal->config.RUN_SEED6;
			seed1[6] = SDglobal->config.RUN_SEED7;
			seed1[7] = SDglobal->config.RUN_SEED8;
			seed1[8] = SDglobal->config.RUN_SEED9;
			seed1[9] = SDglobal->config.RUN_SEED10;
			seed1[10] = SDglobal->config.RUN_SEED11;
			seed1[11] = SDglobal->config.RUN_SEED12;
			seed1[12] = SDglobal->config.RUN_SEED13;
			seed1[13] = SDglobal->config.RUN_SEED14;
			seed1[14] = SDglobal->config.RUN_SEED15;
			seed1[15] = SDglobal->config.RUN_SEED16;
			seed1[16] = SDglobal->config.RUN_SEED17;
			seed1[17] = SDglobal->config.RUN_SEED18;
			seed1[18] = SDglobal->config.RUN_SEED19;
			seed1[19] = SDglobal->config.RUN_SEED20;
			seed1[20] = SDglobal->config.RUN_SEED21;
			seed1[21] = SDglobal->config.RUN_SEED22;
			seed1[22] = SDglobal->config.RUN_SEED23;
			seed1[23] = SDglobal->config.RUN_SEED24;
			seed1[24] = SDglobal->config.RUN_SEED25;
			seed1[25] = SDglobal->config.RUN_SEED26;
			seed1[26] = SDglobal->config.RUN_SEED27;
			seed1[27] = SDglobal->config.RUN_SEED28;
			seed1[28] = SDglobal->config.RUN_SEED29;
			seed1[29] = SDglobal->config.RUN_SEED30;
		}
		else {
			generate_seed(seed1, seed2);
		}
	}
	else {
		if (SDglobal->config.AUTO_SEED == 0) {
			seed1[0] = SDglobal->config.RUN_SEED1;
		}
		else {
			generate_seed(seed1, seed2);
		}
	}

	seed2[0] = SDglobal->config.EVAL_SEED1;

	printf("\nBeginning SD...\n\n");
	/* Solve all of the problems */
	for (n = 0; n < num_probs; n++) {
		/* Find the name of the current problem */
		if (argc > 2 && argc < 5)
			filename_number(fname, 4, 1000, n + start);
		if (load_core_cpx(&probptr, &ident, fname, objsen)) {
			if (load_stoch(SDglobal, probptr, ident, fname)) {
				if (load_time(&row, &col, ident, fname)) {
					sort_omegas(SDglobal, col);
					cipher_omegas(SDglobal);

					if (SDglobal->config.MIN_ITER < SDglobal->config.ITER_FACT * col)
						SDglobal->config.MIN_ITER = SDglobal->config.ITER_FACT * col;

					/* Solve the original combined problem */
					change_solver_primal(probptr); //added by Yifan to change solver to primal simplex
					setup_problem(probptr);
					write_prob(probptr, "orig.lp");
					solve_problem(SDglobal, probptr);

					/* Get the solution from the mean problem */
					if (!(x_k = arr_alloc(probptr->mac+1, double)))
						err_msg("Allocation", "main", "x_k");
					/* Allocate memory for the copy of mean value solution 04/25/2013 Yifan */
					if (!(original_x_k = arr_alloc(probptr->mac+1, double)))
						err_msg("Allocation", "main", "x_k");

					get_primal(x_k, probptr, probptr->mac);
					/* Make a copy of the mean value solution Yifan 04/25/2013 */
					copy_arr(original_x_k, x_k, probptr->mac);

					objective = get_objective(probptr);

					/* First stage solutions to the mean value problem Yifan*/
					print_vect(x_k, probptr->mac, "x_k from the mean problem");

					printf("Objective from the mean problem : %f\n", objective);
					get_lower_bound(SDglobal, probptr, row, col);
					remove_problem(probptr);

					/* Now split the problem and solve it with SD */
					num_rv = SDglobal->omegas.num_omega;
					num_cipher = SDglobal->omegas.num_cipher;
					printf("\nnum_rv = %d;\n", num_rv);
					printf("num_cipher = %d;\n\n", num_cipher);

					if (SDglobal->config.MULTIPLE_REP == 1) {
						for (idx = 1; idx < 2; idx++)
						{/* this layer control the tolerance of SD runs */
							for (cnt = 0; cnt < BATCH_SIZE; cnt++)
							{ /* this layer control the seed used in each SD replication */
								SDglobal->config.EVAL_SEED1 = seed2[0];
								SDglobal->config.RUN_SEED = seed1[cnt];
								SDglobal->store_flag = FALSE;
								SDglobal->pi_flag[0] = FALSE;
								SDglobal->pi_flag[1] = FALSE;
								SDglobal->pi_flag[2] = FALSE;
								/* Take the mean value solution as the initial candidate solution 04/25/2013 Yifan */
								copy_arr(x_k, original_x_k, probptr->mac);
								solve_SD(SDglobal, probptr, x_k, num_rv, num_cipher, row, col, fname, cnt);
								if (cnt == 2 && (SDglobal->average_flag == 1 || SDglobal->obj_flag == 1 ) && SDglobal->config.OVERRIDE == 1) {
									break;
								}

							}
						}
					}
					else {
						SDglobal->config.RUN_SEED = seed1[0];
						solve_SD(SDglobal, probptr, x_k, num_rv, num_cipher, row, col, fname, cnt);
					}

					mem_free(x_k);
					/* Release the copy of mean value solution Yifan 04/25/2013 */
					mem_free(original_x_k);
				} /* End of load_time. */
				else {
					printf("File error in load_time, problem #%d", n);
					return 1;
				}
			} /* End of load_stoch. */
			else {
				printf("File error in load_stoch, problem #%d", n);
				return 1;
			}
		} /* End of load_core_cpx.  zl trial. */
		else {
			printf("Failed in load_core_cpx, problem #%d", n);
			return 1;
		}

		/* Clean up */
		free_ident(ident);
		free_omegas(SDglobal);
	}

	printf("\nEnding SD...\n\n");

	log_stop(SDglobal);
	free_one_prob(probptr);
	freeSDglobal(SDglobal);

	/* Release the CPLEX environment. */
	close_Solver();
#ifdef SD_win
	sd_mv_output_files(fname);
#else
	if (RESUME_FLAG) {
		sd_mv_resume_files(buffer3, buffer2, fname);
	}
	sd_mv_output_files(buffer1, buffer2, fname);
#endif
	return 0;
}

sdglobal_type *setupSDglobal() {
	sdglobal_type *SDglobal;

	SDglobal = (sdglobal_type *) mem_malloc(sizeof(sdglobal_type)); /* All global variables */
	SDglobal->batch_problem = NULL;
	SDglobal->batch_incumb = NULL;
	SDglobal->bcuts = NULL;
	SDglobal->bfcuts = NULL;
	SDglobal->bfcuts_pool = NULL;
	SDglobal->Obj_lb = NULL;
	SDglobal->quad_v = NULL;
	SDglobal->Bbar = NULL;

	/* scalars */
	SDglobal->MEM_USED = 0;
	SDglobal->Abar = 0;
	SDglobal->MALLOC = 0;
	SDglobal->average_flag = 0;
	SDglobal->obj_flag = 0;
	SDglobal->resume_flag = FALSE;
	SDglobal->pi_flag[0] = SDglobal->pi_flag[1] = SDglobal->pi_flag[2] = FALSE;

	return SDglobal;
}//END setupSDglobal()

void freeSDglobal(sdglobal_type *SDglobal) {
	int n;

	if (SDglobal) {
		if (SDglobal->batch_problem) free_one_prob(SDglobal->batch_problem);
		if (SDglobal->Bbar) mem_free(SDglobal->Bbar);
		if (SDglobal->Obj_lb) mem_free(SDglobal->Obj_lb);
		if (SDglobal->quad_v) mem_free(SDglobal->quad_v);
		if (SDglobal->batch_incumb) {
			if ( SDglobal->batch_incumb->incumb_x) {
				for ( n = 0; n < BATCH_SIZE; n++ )
					mem_free(SDglobal->batch_incumb->incumb_x[n]);
				mem_free(SDglobal->batch_incumb->incumb_x);
			}
			if ( SDglobal->batch_incumb->R_Master_pi) {
				for ( n = 0; n < BATCH_SIZE; n++ )
					mem_free(SDglobal->batch_incumb->R_Master_pi[n]);
				mem_free(SDglobal->batch_incumb->R_Master_pi);
			}
			if ( SDglobal->batch_incumb->R_Master_dj) {
				for ( n = 0; n < BATCH_SIZE; n++ )
					mem_free(SDglobal->batch_incumb->R_Master_dj[n]);
				mem_free(SDglobal->batch_incumb->R_Master_dj);
			}
			mem_free(SDglobal->batch_incumb);
		}
		if (SDglobal->bcuts) free_bcuts(SDglobal->bcuts);
		if( SDglobal->bfcuts) free_bcuts(SDglobal->bfcuts);
		if (SDglobal->bfcuts_pool) free_bcuts(SDglobal->bfcuts_pool);
		mem_free(SDglobal);
	}

}

void sd_create_output_folder(sdglobal_type* sd_global, char *buffer1, char *buffer2, char *fname) {
	int status;
	strcpy(buffer1, "mkdir ./sdoutput");
	status = system(buffer1);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcat(buffer1, "/");
	strcat(buffer1, fname);
	status = system(buffer1);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcat(buffer1, "/");
	if (sd_global->config.EPSILON==0.01) {
		strcat(buffer1, "loose");
	}
	else if (sd_global->config.EPSILON==0.001){
		strcat(buffer1, "nominal");
	}
	else{
		strcat(buffer1, "tight");
	}
	status = system(buffer1);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcpy(buffer1, "./sdoutput/");
	strcat(buffer1, fname);
	strcat(buffer1, "/");
	if (sd_global->config.EPSILON==0.01) {
		strcat(buffer1, "loose");
	}
	else if (sd_global->config.EPSILON==0.001){
		strcat(buffer1, "nominal");
	}
	else{
		strcat(buffer1, "tight");
	}

}

void sd_create_resume_folder(sdglobal_type* sd_global, char *buffer3, char *buffer2, char *fname) {
	int status;
	strcpy(buffer3, "mkdir ./sdresume");
	status = system(buffer3);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcat(buffer3, "/");
	strcat(buffer3, fname);
	status = system(buffer3);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcat(buffer3, "/");
	if (sd_global->config.EPSILON==0.01) {
		strcat(buffer3, "loose");
	}
	else if (sd_global->config.EPSILON==0.001){
		strcat(buffer3, "nominal");
	}
	else{
		strcat(buffer3, "tight");
	}
	status = system(buffer3);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
	strcpy(buffer3, "./sdresume/");
	strcat(buffer3, fname);
	strcat(buffer3, "/");
	if (sd_global->config.EPSILON==0.01) {
		strcat(buffer3, "loose");
	}
	else if (sd_global->config.EPSILON==0.001){
		strcat(buffer3, "nominal");
	}
	else{
		strcat(buffer3, "tight");
	}
}

BOOL sd_check_resume_folder(sdglobal_type* sd_global, char *fname) {
	char buffer1[128],buffer2[128];
	int usr_response=-1;
	BOOL tight_flag = FALSE, nominal_flag = FALSE, loose_flag = FALSE;

	strcpy(buffer1, "./sdresume/");
	strcat(buffer1, fname);
	strcat(buffer1, "/");

	/* First check tight tolerance */
	strcpy(buffer2, buffer1);
	tight_flag = sd_check_resume_data(buffer2, "tight");

	if (tight_flag) {
		printf("Data of tight tolerance is available for SD to resume from, type 1 to resume or 0 to start SD from scratch:");
		if (scanf("%d", &usr_response) != 1) {
			printf("Failed to read user_response");
		}
		while ((usr_response!=0)&&(usr_response!=1)) {
			printf("Please type 1 to resume SD or 0 to start SD from scratch:");
			if (scanf("%d", &usr_response) != 1) {
				printf("Failed to read usr_response");
			}
		}
	}

	/* Then check nominal tolerance */
	strcpy(buffer2, buffer1);
	nominal_flag = sd_check_resume_data(buffer2, "nominal");

	if (nominal_flag) {
		printf("Data of nominal tolerance is available for SD to resume, type 1 to resume or 0 to start from scratch:");
		if (scanf("%d", &usr_response) != 1) {
			printf("Failed to read user_response");
		}
		while ((usr_response!=0)&&(usr_response!=1)) {
			printf("Please type 1 to resume SD or 0 to start SD from scratch:");
			if (scanf("%d", &usr_response) != 1) {
				printf("Failed to read usr_response");
			}
		}
	}

	/* Lastly, check loose tolerance */
	strcpy(buffer2, buffer1);
	loose_flag = sd_check_resume_data(buffer2, "loose");

	if (loose_flag) {
		printf("Data of loose tolerance is available for SD to resume, type 1 to resume or 0 to start from scratch:");
		if (scanf("%d", &usr_response) != 1) {
			printf("Failed to read user_response");
		}
		while ((usr_response!=0)&&(usr_response!=1)) {
			printf("Please type 1 to resume SD or 0 to start SD from scratch:");
			if (scanf("%d", &usr_response) != 1) {
				printf("Failed to read usr_response");
			}
		}
	}

	if (usr_response == 1) {
		return TRUE;
	}
	else
	{
		return FALSE;
	}

}

BOOL sd_check_resume_data(char *buffer2, char *tolerance) {
	BOOL flag=FALSE;
	char buffer3[128];
	char temp_buffer[128];
	int cnt;
	char rep_number[16];

	strcat(buffer2, tolerance);
	strcat(buffer2, "/");
	strcpy(buffer3, buffer2);
	strcat(buffer2, "resume_data");
	strcat(buffer3, "resume");

	for(cnt=0; cnt < BATCH_SIZE; cnt++) {
		sprintf(rep_number, "%d", cnt);

		strcpy(temp_buffer, buffer2);
		strcat(temp_buffer, rep_number);
		strcat(temp_buffer, ".txt");

		if( access( temp_buffer, F_OK ) != -1 ) {
			// file exists
		} else {
			break;
		}

		strcpy(temp_buffer, buffer3);
		strcat(temp_buffer, rep_number);
		strcat(temp_buffer, ".lp");

		if( access( temp_buffer, F_OK ) != -1 ) {
			// file exists
		} else {
			break;
		}

		if (cnt==(BATCH_SIZE-1)) {
			flag = TRUE;
		}
	}

	return flag;
}


#ifdef SD_win
void sd_mv_output_files(char *fname)
{
	TCHAR path[BUFFER_SIZE];
	TCHAR inst_dir[BUFFER_SIZE];
	TCHAR buff[BUFFER_SIZE];
	char file_name[BUFFER_SIZE];
	TCHAR L_file_name[BUFFER_SIZE];

	if(!GetCurrentDirectory(BUFFER_SIZE, path))
		printf("GetCurrentDirectory() failed!\n");

	printf("Your current directory is: %S\n", path);

	wcscat(path,L"\\sdoutput");
	printf("Your SD Output directory is: %S\n", path);
	CreateDirectory(path,NULL);

	mbstowcs(inst_dir,fname,NAME_SIZE);
	printf("Your instance's name is: %S\n", inst_dir);
	wcscat(path,L"\\");
	wcscat(path,inst_dir);
	printf("Your instance SD Result directory is: %S\n", path);
	CreateDirectory(path,NULL);

	wcscpy(buff, path);
	wcscat(buff, L"\\after_coef_change.lp");
	MoveFile(L"after_coef_change.lp", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_dual.out");
	MoveFile(L"Batch_dual.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_Obj.out");
	MoveFile(L"Batch_Obj.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_x.out");
	MoveFile(L"Batch_x.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\before_coef_change.lp");
	MoveFile(L"before_coef_change.lp", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\eval.dat");
	MoveFile(L"eval.dat", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\final_master.lp");
	MoveFile(L"final_master.lp", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\final-batch-prob.lp");
	MoveFile(L"final-batch-prob.lp", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\incumb.out");
	MoveFile(L"incumb.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Master_Dual.out");
	MoveFile(L"Master_Dual.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Master_Obj.out");
	MoveFile(L"Master_Obj.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\orig.lp");
	MoveFile(L"orig.lp", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\Summary.out");
	MoveFile(L"Summary.out", buff);

	wcscpy(buff, path);
	wcscat(buff, L"\\time_sample.out");
	MoveFile(L"time_sample.out", buff);

	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".detailed_rep_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".detailed_rep_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);

	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".detailed_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".detailed_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);

	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".obj.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".obj.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);

	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".time.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".time.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);
}
#else
void sd_mv_output_files(char *buffer1, char *buffer2, char *fname) {
	int status;
	if (RESUME_FLAG) {
		strcpy(buffer2, "mv resume_time store_time *.out *.lp *.dat ");
	}
	else{
		strcpy(buffer2, "mv *.out *.lp *.dat ");
	}

	strcat(buffer2, buffer1);
	status = system(buffer2);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
}

void sd_mv_resume_files(char *buffer3, char *buffer2, char *fname) {
	int status;
	strcpy(buffer2, "mv resume_data*.txt resume*.lp ");
	strcat(buffer2, buffer3);
	status = system(buffer2);
	if(status == -1){
		printf("system() call fails.\n");
		exit(1);
	}
}
#endif

//mem_free(sd_global->batch_incumb->incumb_x);
//mem_free(sd_global->batch_incumb->R_Master_pi);
//mem_free(sd_global->batch_incumb->R_Master_dj);
//mem_free(sd_global->batch_incumb);
//mem_free(sd_global->Obj_lb);
//mem_free(sd_global->quad_v);
//
//free_bcuts(sd_global->bcuts);
//free_bcuts(sd_global->bfcuts_pool);
//for (idx=0; idx<BATCH_SIZE; idx++) {
//	if (sd_global->bfcuts->batch[idx]!= NULL) {
//		if (sd_global->bfcuts->batch[idx]->val!=NULL) {
//			mem_free(sd_global->bfcuts->batch[idx]->val);
//		}
//	}
//	if (sd_global->bfcuts->batch[idx]!=NULL) {
//		mem_free(sd_global->bfcuts->batch[idx]);
//	}
//}
//
//mem_free(sd_global->bfcuts->batch);
//mem_free(sd_global->bfcuts);
//
//free_one_prob(sd_global->batch_problem);

