/************************************************************************
 * bdmatch.h
 *
 * Main header file for bdmatch bp code.
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
 ************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <string.h>

#define MAX_ITER 1000
#define INTERVAL 2
#define NEG_INF -999999
#define CIRC_BUFF 100
#define CONVERG_THR 1e-4
#define INF 1e16

#ifndef __UPDATEARGS
#define __UPDATEARGS
typedef struct {
    float **oldB, **B, **W;
    int *UB, *LB, **backinds, **inds, *rowcounts,N;
    int *myRows, myCount;
} updateArgs;
#endif


void *updateB(void *args);

int quickselect(float *V, int N, int k);

int recursiveSelect(float *V, int *inds, int start, int end, int k);

extern float damping;
extern int THREADS;
