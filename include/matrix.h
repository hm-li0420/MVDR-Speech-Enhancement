#ifndef COMPLEX_MATRIX
#define COMPLEX_MATRIX
#include "kiss_fftr.h"

typedef struct Complex_Matrix_State {
	// ������к���
	int row, column; 
	kiss_fft_cpx* arrayComplex;
}complex_matrix_state;

typedef struct Matrix_State {
	// ������к���
	int row, column;
	float* array;
}matrix_state;

// ��ʼ������
void InitComplex(kiss_fft_cpx* Complex);

// ��ʼ���������� row������ column�� ����
void InitComplexMatrix(complex_matrix_state* matrix, int row, int column);

// ��ʼ���������� row������ column�� ����
void InitMatrix(matrix_state* matrix, int row, int column);

// �����������
complex_matrix_state* cpx_matmul(complex_matrix_state* cpx_matrix_a, complex_matrix_state* cpx_matrix_b);

// �����������
kiss_fft_cpx* cpx_divide(kiss_fft_cpx* a, kiss_fft_cpx* b);

// �����������
matrix_state* matmul(matrix_state* matrix_a, matrix_state* matrix_b);

// ����������ת��
complex_matrix_state* cpx_matrix_transpose(complex_matrix_state* cpx_matrix);

// ������������
complex_matrix_state* cpx_matrix_inverse(complex_matrix_state* cpx_matrix);

// �ͷų�������
void matrix_destroy(matrix_state* matrix);

// �ͷŸ�������
void cpx_matrix_destroy(complex_matrix_state* cpx_matrix);
#endif