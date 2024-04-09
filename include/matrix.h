#ifndef COMPLEX_MATRIX
#define COMPLEX_MATRIX
#include "kiss_fftr.h"

typedef struct Complex_Matrix_State {
	// 矩阵的行和列
	int row, column; 
	kiss_fft_cpx* arrayComplex;
}complex_matrix_state;

typedef struct Matrix_State {
	// 矩阵的行和列
	int row, column;
	float* array;
}matrix_state;

// 初始化复数
void InitComplex(kiss_fft_cpx* Complex);

// 初始化复数矩阵 row：行数 column： 列数
void InitComplexMatrix(complex_matrix_state* matrix, int row, int column);

// 初始化常数矩阵 row：行数 column： 列数
void InitMatrix(matrix_state* matrix, int row, int column);

// 复数矩阵相乘
complex_matrix_state* cpx_matmul(complex_matrix_state* cpx_matrix_a, complex_matrix_state* cpx_matrix_b);

// 复数矩阵相除
kiss_fft_cpx* cpx_divide(kiss_fft_cpx* a, kiss_fft_cpx* b);

// 常数矩阵相乘
matrix_state* matmul(matrix_state* matrix_a, matrix_state* matrix_b);

// 复数矩阵共轭转置
complex_matrix_state* cpx_matrix_transpose(complex_matrix_state* cpx_matrix);

// 复数矩阵求逆
complex_matrix_state* cpx_matrix_inverse(complex_matrix_state* cpx_matrix);

// 释放常数矩阵
void matrix_destroy(matrix_state* matrix);

// 释放复数矩阵
void cpx_matrix_destroy(complex_matrix_state* cpx_matrix);
#endif