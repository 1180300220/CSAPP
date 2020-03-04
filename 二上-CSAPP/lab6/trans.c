/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i,j;
    if(M == 32 && N == 32){
        for(i = 0;i < M;i = i+8){
            for(j = 0;j < N;j++){
                int a1 = A[j][i];
                int a2 = A[j][i+1];
                int a3 = A[j][i+2];
                int a4 = A[j][i+3];
                int a5 = A[j][i+4];
                int a6 = A[j][i+5];
                int a7 = A[j][i+6];
                int a8 = A[j][i+7];
                B[i][j] = a1;
                B[i+1][j] = a2;
                B[i+2][j] = a3;
                B[i+3][j] = a4;
                B[i+4][j] = a5;
                B[i+5][j] = a6;
                B[i+6][j] = a7;
                B[i+7][j] = a8;
            }
        }
    }
    if(M == 64 && N == 64){
        for(int i = 0;i < N;i += 8){
            for(int j = 0;j < M;j += 8){
                for(int k = i;k < i+4;k++){
                    //首先读取1，2，并将其放在B中的1，3的位置
                    int a1 = A[k][j];
                    int a2 = A[k][j+1];
                    int a3 = A[k][j+2];
                    int a4 = A[k][j+3];
                    int a5 = A[k][j+4];
                    int a6 = A[k][j+5];
                    int a7 = A[k][j+6];
                    int a8 = A[k][j+7];
                    B[j][k] = a1;
                    B[j+1][k] = a2;
                    B[j+2][k] = a3;
                    B[j+3][k] = a4;
                    B[j][k+4] = a8;
                    B[j+1][k+4] = a7;
                    B[j+2][k+4] = a6;
                    B[j+3][k+4] = a5;
                }
                for(int p = 0;p < 4; p++){
                    
                    //首先将3，4取出
                    int a1 = A[i+4][j+3-p];
                    int a2 = A[i+5][j+3-p];
                    int a3 = A[i+6][j+3-p];
                    int a4 = A[i+7][j+3-p];
                    int a5 = A[i+4][j+4+p];
                    int a6 = A[i+5][j+4+p];
                    int a7 = A[i+6][j+4+p];
                    int a8 = A[i+7][j+4+p];
                    //然后将A中的2（此时在B中的3）移动到B中的2.
                    B[j+4+p][i] = B[j+3-p][i+4];
                    B[j+4+p][i+1] = B[j+3-p][i+5];
                    B[j+4+p][i+2] = B[j+3-p][i+6];
                    B[j+4+p][i+3] = B[j+3-p][i+7];
                    //将A中3，4（刚才取出的）放入对应位置
                    B[j+3-p][i+4] = a1;
                    B[j+3-p][i+5] = a2;
                    B[j+3-p][i+6] = a3;
                    B[j+3-p][i+7] = a4;
                    B[j+4+p][i+4] = a5;
                    B[j+4+p][i+5] = a6;
                    B[j+4+p][i+6] = a7;
                    B[j+4+p][i+7] = a8;
                    
                }
            }
        }
    }
    if(M == 61 && N == 67){
        for(int i = 0;i < N;i+=16){
            for(int j = 0;j < M;j+= 16){
                for(int k = i;k < i+16 && k < N;++k){
                    int if_Diagonal = -1;
                    int Diagonal_tmp = 0;
                    for(int p = j;p <j + 16 && p < M;++p){
                        //如果是对角线上的元素，则将其单独提出。
                        if(p == k){
                            Diagonal_tmp = A[p][p];
                            if_Diagonal = p;
                        //如果不是
                        }else{
                            B[p][k] = A[k][p];
                        }
                    }
                    //将对角线上的元素单独放入
                    if(if_Diagonal != -1){
                        B[if_Diagonal][if_Diagonal] = Diagonal_tmp;
                    }
                }
            }
        }
    }
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

