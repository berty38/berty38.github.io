/************************************************************************
 * bp.c
 *
 * belief update function. 
 *
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

float damping;

void *updateB(void *args)
{
    int c,i,j,k,poscount,bth,bplus,dth,dplus;
    updateArgs *arg = args;
    float **oldB = arg->oldB, **B = arg->B, **W = arg->W;
    int *UB = arg->UB, *LB = arg->LB, **backinds = arg->backinds,
    **inds = arg->inds, *rowcounts = arg->rowcounts,
    *myRows = arg->myRows, myCount = arg->myCount;
    
    for (c=0; c<myCount; c++) {
        
        j = myRows[c];
        
        poscount = 0;
        for (i=0; i<rowcounts[j]; i++) {
            
            if (oldB[j][i]>0)
                poscount++;
        }
        
        if (UB[j]==0) {
            for (i=0; i<rowcounts[j]; i++) {
                k = inds[j][i];
                B[k][backinds[j][i]] = -INF;
            }
        } else if ((LB[j]<poscount && poscount<UB[j]) ||
        (LB[j]==0 && poscount==0)) {
            for (i=0; i<rowcounts[j]; i++) {
                k = inds[j][i];
                B[k][backinds[j][i]] = W[k][backinds[j][i]]+W[j][i];
            }
        } else {
            
            bth = quickselect(oldB[j],rowcounts[j],LB[j]-1);
            bplus = quickselect(oldB[j],rowcounts[j],LB[j]);
            dth = quickselect(oldB[j],rowcounts[j],UB[j]-1);
            dplus = quickselect(oldB[j],rowcounts[j],UB[j]);
            
            for (i=0; i<rowcounts[j]; i++) {
                k = inds[j][i];
                
                if (poscount<=LB[j]) {
                    if (oldB[j][i]>=oldB[j][bth]) {
                        if (bplus<0) {
                            B[k][backinds[j][i]] = INF;
                        } else {
                            B[k][backinds[j][i]] =
                            W[k][backinds[j][i]]+W[j][i]-oldB[j][bplus];
                        }
                    } else {
                        if (poscount==LB[j] && UB[j]>LB[j]) {
                            B[k][backinds[j][i]] =
                            W[k][backinds[j][i]]+W[j][i];
                        } else {
                            B[k][backinds[j][i]] =
                            W[k][backinds[j][i]]+W[j][i]-oldB[j][bth];
                        }
                    }
                } else if (poscount==UB[j]) {
                    if (oldB[j][i]>=oldB[j][dth]) {
                        B[k][backinds[j][i]] =
                        W[k][backinds[j][i]]+W[j][i];
                    } else {
                        B[k][backinds[j][i]] =
                        W[k][backinds[j][i]]+W[j][i]-oldB[j][dth];
                    }
                } else if (poscount > UB[j]) {
                    if (oldB[j][i]>=oldB[j][dth]) {
                        if (dplus<0)
                            B[k][backinds[j][i]] = INF;
                        else
                            B[k][backinds[j][i]] =
                            W[k][backinds[j][i]]+W[j][i]-oldB[j][dplus];
                    } else {
                        B[k][backinds[j][i]] =
                        W[k][backinds[j][i]]+W[j][i]-oldB[j][dth];
                    }
                }
                
                /*
                 if (isnan(B[k][backinds[j][i]])) {
                    printf("NaN detected: B[%d][%d]",i,backinds[j][i]);
                    printf("poscount = %d\n",poscount);
                    printf("LB = %d, UB = %d\n", LB[j], UB[j]);
                }
                 */
                
                B[k][backinds[j][i]] = damping*B[k][backinds[j][i]] 
                + (1-damping)*oldB[k][backinds[j][i]];
            }
        }
    }
    free(args);
    
    return NULL;
}
