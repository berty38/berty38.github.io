/************************************************************************
 * bdmatch.c
 *
 * Main function for bdmatch bp code.
 *
 * Copyright (C) 2008  Bert Huang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/


#include "bdmatch.h"

int THREADS = 4;

int main(int argc, char *argv[])
{
    int i,j,k,N,*LB,*UB,*colsum,iter=1,converged, poscount,**inds,
    *rowcounts,**backinds,*colcounts,myCount[THREADS], **myRows,
    *updatethread, cbuffpos=0, oscillation = 10, row, col, nnz=0,
    iterchars, maxiter = MAX_ITER, verbose = 1, sparse = 1, flags = 0;
    float **W, **B, **oldB, change, *rowsums, *oldrowsums, val, rowf, colf,
    lbf, ubf, threshold, cbuff[CIRC_BUFF], Nf, nnzf;
    FILE *matrix=NULL, *bounds=NULL,*out=NULL;
    updateArgs *arg;
    pthread_t *beliefUpdates;
    
    /* load arguments from command line */
    damping = 1;
    
    i = 1;
    while(i<argc) {
        if (strcmp(argv[i],"-w")==0) {
            if (argc>i+1) {
                matrix = fopen(argv[i+1],"r");
                i+=2;
                flags++;
            } else {
                printf("Error, a filename must follow '-w'\n");
                return 0;
            }
        } else if (strcmp(argv[i], "-o")==0) {
            if (argc>i+1) {
                out = fopen(argv[i+1],"w");
                i+=2;
                flags++;
            } else {
                printf("Error, a filename must follow '-o'\n");
                return 0;
            }
        } else if (strcmp(argv[i], "-d")==0) {
            if (argc>i+1) {
                bounds = fopen(argv[i+1],"r");
                i+=2;
                flags++;
            } else {
                printf("Error, a filename must follow '-d'\n");
                return 0;
            }
        } else if (strcmp(argv[i], "-i")==0) {
            if (argc>i+1) {
                maxiter = atoi(argv[i+1]);
                i+=2;
            } else {
                printf("Error, a number must follow '-i'\n");
                return 0;
            }
        } else if (strcmp(argv[i], "-damp")==0) {
            if (argc>i+1) {
                damping = atof(argv[i+1]);
                i+=2;
            } else {
                printf("Error, a number must follow '-damp'\n");
            }
        } else if (strcmp(argv[i], "-v")==0) {
            if (argc>i+1) {
                verbose = atoi(argv[i+1]);
                i+=2;
            } else {
                printf("Error, a number must follow '-v'\n");
            }
        } else if (strcmp(argv[i], "-s")==0) {
            if (argc>i+1) {
                sparse = atoi(argv[i+1]);
                i+=2;
            } else {
                printf("Error, a number must follow '-s'\n");
            }
	} else if (strcmp(argv[i], "-t")==0) {
	  if (argc>i+1) {
	    THREADS = atoi(argv[i+1]);
	    i+=2;
	  } else {
	    printf("Error, a number must follow '-t'\n");
	  }
        } else if (strcmp(argv[i], "-h")==0) {
	    printf("bdmatch version 0.1b\n");
            printf("USAGE: bdmatch [-arg1 val] [-arg2 val] [..]\n\n");
            printf("Flags:\n\n");
            printf("\t-w [filename]\tWeight file\n");
            printf("\t-d [filename]\tDegree file\n");
            printf("\t-o [filename]\tOutput file\n");
            printf("\t-i [integer]\tMaximum iterations (default 1000)\n");
            printf("\t-damp [float]\tDamping proportion (default 1)\n");
            printf("\t-h\t\tHelp\n");
            printf("\t-s [0 1]\tSparse input (default 1)\n");
            printf("\t-t [integer]\tNumber of threads (default 4)\n");
            printf("\t-v [0 1 2]\tVerbose (default 1)\n");
            printf("\n\nSparse format should have 'N nnz' as first line\n");
            
            printf("bdmatch  Copyright (C) 2008 Bert Huang\n");
            printf("This program comes with ABSOLUTELY NO WARRANTY. ");
            printf("This is free software, and you are welcome to");
	    printf(" redistribute it under certain conditions. ");
	    printf("See LICENSE for details.\n");
            i++;
        } else {
            printf("Error, unrecognized command line argument\n");
            return 0;
        }
    }
    
    if (flags<3) {
        printf("Error, need at least weights, degrees and output filenames\n");
        return 0;
    }
    
    if (verbose>=1) {
        printf("Verbose level is %d. ", verbose);
        printf("Damping factor is %f. ", damping);
	printf("Running with %d threads.\n", THREADS);
    }
    
    fscanf(matrix, "%f %f\n", &Nf,&nnzf);
    nnz = (int)nnzf;
    N = (int)Nf;
   
    if (verbose>=1)
        printf("%d nodes, %d edges\n", N,nnz);
    
  /* allocate memory */
    rowcounts = (int*)malloc(N*sizeof(int));
    rowsums = (float*)malloc(N*sizeof(float));
    oldrowsums = (float*)malloc(N*sizeof(float));
    colcounts = (int*)malloc(N*sizeof(int));
    backinds = (int**)malloc(N*sizeof(int*));
    W = (float**)malloc(N*sizeof(float*));
    B = (float**)malloc(N*sizeof(float*));
    oldB = (float**)malloc(N*sizeof(float*));
    inds = (int**)malloc(N*sizeof(int*));
    
  /* first count nonzeros */
    for (i=0; i<N; i++)
        rowcounts[i]=0;
    
    for (i=0; i<nnz; i++) {
        fscanf(matrix, "%f %f %f\n", &rowf,&colf,&val);
        row = (int)rowf-1;
        rowcounts[row]++;
    }
    
    rewind(matrix);
    
    for (i=0; i<N; i++) {
    /* printf("rowcounts[%d]=%d\n",i,rowcounts[i]); */
        
        inds[i] = malloc(rowcounts[i]*sizeof(int));
        backinds[i] = malloc(rowcounts[i]*sizeof(int));
        W[i] = malloc(rowcounts[i]*sizeof(float));
        B[i] = malloc(rowcounts[i]*sizeof(float));
        oldB[i] = malloc(rowcounts[i]*sizeof(float));
        
        rowcounts[i]=0; /*reset rowcounts to use for indexing */
    }
    
    for (j=0; j<N; j++) {
        colcounts[j]=0;
    }
    iterchars = 0;
    
    fscanf(matrix, "%f %f\n", &Nf,&nnzf);
    nnz = (int)nnzf;
    N = (int)Nf;
       
    for (i=0; i<nnz; i++) {
        fscanf(matrix, "%f %f %f\n", &rowf,&colf,&val);
        row = (int)rowf-1;
        col = (int)colf-1;
        
        if (verbose>=2)
            printf("Connecting node %d to %d with weight %f\n", row, col, val);
        
        W[row][rowcounts[row]] = val;
        B[row][rowcounts[row]] = val;
        oldB[row][rowcounts[row]] = val;
        inds[row][rowcounts[row]] = col;
        backinds[row][rowcounts[row]] = colcounts[col]++;
        rowcounts[row]++;
    }
    
    fclose(matrix);
    
  /* allocate and read in LB */
    
    LB = (int*)malloc(N*sizeof(int));
    UB = (int*)malloc(N*sizeof(int));
    
    for (i=0; i<N; i++) {
        fscanf(bounds, "%f %f\n", &lbf, &ubf);
        LB[i] = (int)lbf;
        UB[i] = (int)ubf;
        
        if (verbose>=3)
            printf("bounds for %d, lb %d, ub %d\n", i, LB[i], UB[i]);
    }
    fclose(bounds);
    
  /*validate inputs */
    
    for (i=0; i<N; i++) {
        if (LB[i]>UB[i]) {
            printf("Lower bound cannot be greater than upper\n");
            return 0;
        }
    }
    
    if (verbose>=1)
        printf("Data loaded. Starting BP\n");
    
    colsum = (int*)malloc(N*sizeof(int));
    
  /*split up threads */
    
    for (i=0; i<THREADS; i++)
        myCount[i]=0;
    
    updatethread = (int*)malloc(N*sizeof(int));
    
    for (i=0; i<N; i++) {
        updatethread[i] = i%THREADS;
        myCount[updatethread[i]]++;
    }
    
    myRows = (int**)malloc(THREADS*sizeof(int*));
    
    for (i=0; i<THREADS; i++) {
        myRows[i] = (int*)malloc(myCount[i]*sizeof(int));
        k=0;
        for (j=0; j<N; j++) {
            if (updatethread[j]==i)
                myRows[i][k++] = j;
        }
    }
    
    beliefUpdates = (pthread_t*)malloc(THREADS*sizeof(pthread_t));
    
  /* Run Belief Revision */
    converged = 0;
    change = 1;
    if (verbose>=1)
        printf("Iteration ");
    
    iterchars = 0;
    while(converged<N*INTERVAL && change>CONVERG_THR && oscillation>0) {
    /* copy B */
        for (i=0; i<N; i++) {
            for (j=0; j<rowcounts[i]; j++) {
                oldB[i][j] = B[i][j];
            }
        }
        
    /* update B */
        
        for (i=0; i<THREADS; i++) {
            arg = (updateArgs*)malloc(sizeof(updateArgs));
            
            arg->oldB = oldB;
            arg->B = B;
            arg->W = W;
            arg->LB = LB;
            arg->UB = UB;
            arg->rowcounts = rowcounts;
            arg->inds = inds;
            arg->backinds = backinds;
            arg->N = N;
            
            arg->myCount = myCount[i];
            arg->myRows = myRows[i];
            
            pthread_create(&beliefUpdates[i],NULL,updateB,arg);
            
        }
        
        for (i=0; i<THREADS; i++) {
            pthread_join(beliefUpdates[i],NULL);
        }
        
        if (verbose>=1) {
            while(iterchars>0) {
                printf("\b");
                iterchars--;
            }
            printf("%d",iter);
            iterchars = 1+log10(iter);
            fflush(stdout);
        }
    /* check for convergence */
        
        if (iter%INTERVAL==0) {
      /* track changes */
            
            cbuff[cbuffpos] = 0;
            for (i=0; i<N; i++) {
                if (rowcounts[i]>0) {
                    cbuff[cbuffpos] += fabs(B[i][0]);
                }
            }
            
            for (i=0; i<CIRC_BUFF; i++) {
                if (i!=cbuffpos && fabs(cbuff[i]-cbuff[cbuffpos])<CONVERG_THR) {
                    oscillation--;
                    break;
                }
            }
            
            cbuffpos++;
            if (cbuffpos>=CIRC_BUFF)
                cbuffpos = 0;
            
            
            
            for (i=0; i<N; i++) {
                rowsums[i] = 0;
                oldrowsums[i] = 0;
                for (j=0; j<rowcounts[i]; j++) {
                    if (!isinf(exp(B[i][j])))
                        rowsums[i]+=exp(B[i][j]);
                    if (!isinf(exp(oldB[i][j])))
                        oldrowsums[i]+=exp(oldB[i][j]);
                }
            }
            
            change = 0;
            for (i=0; i<N; i++) {
                if (rowsums[i]==0) {
                    rowsums[i]=1;
                }
                if (oldrowsums[i]==0)
                    oldrowsums[i]=1;
                for (j=0; j<rowcounts[i]; j++) {
                    if ((isinf(exp(oldB[i][j])) && !isinf(exp(B[i][j]))) ||
                    (!isinf(exp(oldB[i][j])) && isinf(exp(B[i][j]))))
                        change+=1;
                    else if (isinf(exp(oldB[i][j])) && isinf(exp(B[i][j])))
                        change = change;
                    else
                        change += fabs(exp(oldB[i][j])/oldrowsums[i]-
                        exp(B[i][j])/rowsums[i]);
                }
            }
            if (isnan(change)) {
                printf("change is NaN! BP will quit but solution ");
                printf("could be invalid. Problem may be infeasible.\n");
            }
        }
        
    /* check for maximum iterations */
        
        if (iter+1>maxiter) {
            if (verbose>=1)
                printf("\nReached Maximum iterations.\n");
            converged = 9999999;
        }
        
        iter++;
        
    }
    
    if (verbose>=1) {
        printf("\n");
        if (change<=CONVERG_THR) {
            printf("Converged to stable beliefs in %d iterations\n",iter);
        } else if (oscillation<1) {
            printf("Stopped after reaching oscillation. \n");
            printf("No feasible solution found or there are multiple maxima. ");
            printf("Outputting best approximate solution. Try damping.\n");
        }
    }
    
    /* output P */
    /*output P and B as sparse matrices */
        
    for (i=0; i<N; i++) {
        
        poscount = 0;
        for (j=0; j<rowcounts[i]; j++)
            if (B[i][j]>0)
                poscount++;
        if (LB[i]>=rowcounts[i]) {
            threshold = -INF;
        } else if (UB[i]<1) {
            threshold = INF;
        } else if (poscount>=UB[i]) {
            threshold = B[i][quickselect(B[i],rowcounts[i],UB[i]-1)];
        } else if (poscount<=LB[i]) {
            if (LB[i]==0)
                threshold = INF;
            else
                threshold = B[i][quickselect(B[i],rowcounts[i],LB[i]-1)];
        } else {
            threshold = B[i][quickselect(B[i],rowcounts[i],poscount-1)];
        }
        
        colcounts[i] = 0;
        
        for (j=0; j<rowcounts[i]; j++) {
            fprintf(out,"%d %d %d %f\n", inds[i][j],i,(B[i][j]>=threshold),B[i][j]);
            if (verbose>=2)
                printf("%d %d %d %f\n", inds[i][j],i,(B[i][j]>=threshold),B[i][j]);
        }
    }
    fclose(out);
    
    /* clean up */
    for (i=0; i<N; i++) {
        free(W[i]);
        free(B[i]);
        free(oldB[i]);
        free(inds[i]);
        free(backinds[i]);
    }
    free(W);
    free(B);
    free(inds);
    free(backinds);
    free(oldB);
    free(LB);
    free(colcounts);
    free(colsum);
    free(rowsums);
    free(oldrowsums);
    free(updatethread);
    return 0;
}

