/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <mex.h>

void jisotonic5(int N, double* BB, double* MSE, double* AA, double* WW)
{
    if (N < 1)
        return;

    double* unweightedcount = (double*)mxMalloc(sizeof(double) * N);
    double* count = (double*)mxMalloc(sizeof(double) * N);
    double* sum = (double*)mxMalloc(sizeof(double) * N);
    double* sumsqr = (double*)mxMalloc(sizeof(double) * N);
    int last_index = -1;

    last_index++;
    unweightedcount[last_index] = 1;
    count[last_index] = WW[0];
    sum[last_index] = AA[0] * WW[0];
    sumsqr[last_index] = AA[0] * AA[0] * WW[0];
    MSE[0] = 0;

    for (int j = 1; j < N; j++) {
        last_index++;
        unweightedcount[last_index] = 1;
        count[last_index] = WW[j];
        sum[last_index] = AA[j] * WW[j];
        sumsqr[last_index] = AA[j] * AA[j] * WW[j];
        MSE[j] = MSE[j - 1];
        while (true) {
            if (last_index <= 0)
                break;
            if (sum[last_index - 1] / count[last_index - 1] < sum[last_index] / count[last_index]) {
                break;
            }
            else {
                double prevMSE = sumsqr[last_index - 1] - sum[last_index - 1] * sum[last_index - 1] / count[last_index - 1];
                prevMSE += sumsqr[last_index] - sum[last_index] * sum[last_index] / count[last_index];
                unweightedcount[last_index - 1] += unweightedcount[last_index];
                count[last_index - 1] += count[last_index];
                sum[last_index - 1] += sum[last_index];
                sumsqr[last_index - 1] += sumsqr[last_index];
                double newMSE = sumsqr[last_index - 1] - sum[last_index - 1] * sum[last_index - 1] / count[last_index - 1];
                MSE[j] += newMSE - prevMSE;
                last_index--;
            }
        }
    }

    int ii = 0;
    for (int k = 0; k <= last_index; k++) {
        for (int cc = 0; cc < unweightedcount[k]; cc++) {
            BB[ii + cc] = sum[k] / count[k];
        }
        ii += unweightedcount[k];
    }
}

// [B,MSE]=jisotonic5_mex(A,weights)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    if (nrhs < 2) {
        mexPrintf("Not enough inputs.\n");
        return;
    }
    if (nlhs < 2) {
        mexPrintf("Not enough outputs.\n");
        return;
    }
    int M = mxGetM(prhs[0]);
    int N = mxGetN(prhs[0]);
    if (M != 1) {
        mexPrintf("Input must be a row vector.\n");
        return;
    }
    if (mxGetN(prhs[1]) != N) {
        mexPrintf("Inconsistent dimensions between A and W.\n");
        return;
    }

    double* AA = mxGetPr(prhs[0]);
    double* WW = mxGetPr(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);

    double* BB = mxGetPr(plhs[0]);
    double* MSE = mxGetPr(plhs[1]);

    jisotonic5(N, BB, MSE, AA, WW);
}
