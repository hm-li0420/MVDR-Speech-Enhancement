#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "bf.h"

void InitComplex(kiss_fft_cpx* Complex)
{
	Complex->r = 0.0;
	Complex->i = 0.0;
}

void InitComplexMatrix(complex_matrix_state* matrix, int row, int column)
{
	int size = row * column * sizeof(kiss_fft_cpx);
	//int size = row * column;
	if (size <= 0)
	{
		puts("ERROE: An invalid matrix!\n");
		return;
	}

	//matrix->arrayComplex = malloc(size);
	matrix->arrayComplex = matrix + sizeof(matrix);
	if (size)
	{
		matrix->row = row;                           			 
		matrix->column = column;
		memset(matrix->arrayComplex, 0, size);
		//for (int index = 0; index < row * column; index++)       
		//{
		//	(matrix->arrayComplex + index)->r = 0.0;
		//	(matrix->arrayComplex + index)->i = 0.0;
		//	//InitComplex(matrix->arrayComplex + index);
		//}
	}
}

void InitMatrix(matrix_state* matrix, int row, int column)
{
	int size = row * column * sizeof(float);
	if (size <= 0)
	{
		puts("ERROE: An invalid matrix!\n");
		return;
	}
	matrix->array = malloc(size);

	if (matrix->array)
	{
		matrix->row = row;
		matrix->column = column;
		memset(matrix->array, 0, size);
	}
}

kiss_fft_cpx* cpx_divide(kiss_fft_cpx* a, kiss_fft_cpx* b)
{
	float norm = b->r * b->r + b->i * b->i;
	kiss_fft_cpx* c = malloc(sizeof(kiss_fft_cpx));
	kiss_fft_cpx* cpx_divide(kiss_fft_cpx * a, kiss_fft_cpx * b);
	
	c->r = (a->r * b->r + a->i * b->i) / norm;
	c->i = (a->i * b->r - a->r * b->i) / norm;
	return c;
}

complex_matrix_state *cpx_matmul(complex_matrix_state* cpx_matrix_a, complex_matrix_state* cpx_matrix_b)
{
	complex_matrix_state* cpx_matrix_c = (complex_matrix_state*)malloc(sizeof(complex_matrix_state));
	int row_i, column_j, ij;
	int idx_a, idx_b, idx_c;
	kiss_fft_cpx temp_cpx;
	temp_cpx.r = 0.f;
	temp_cpx.i = 0.f;
	if (cpx_matrix_a->column != cpx_matrix_b->row)
	{
		puts("ERROE: An incompatable matrix!\n");
		return;
	}
	else
	{
		InitComplexMatrix(cpx_matrix_c, cpx_matrix_a->row, cpx_matrix_b->column);
		for (row_i = 0; row_i < cpx_matrix_c->row; row_i++)
		{
			for (column_j = 0; column_j < cpx_matrix_c->column; column_j++)
			{
				for (ij = 0; ij < cpx_matrix_a->column; ij++)
				{
					idx_a = row_i * cpx_matrix_a->column + ij;
					idx_b = ij * cpx_matrix_b->column + column_j;
					temp_cpx.r += cpx_matrix_a->arrayComplex[idx_a].r * cpx_matrix_b->arrayComplex[idx_b].r - cpx_matrix_a->arrayComplex[idx_a].i * cpx_matrix_b->arrayComplex[idx_b].i;
					temp_cpx.i += cpx_matrix_a->arrayComplex[idx_a].r * cpx_matrix_b->arrayComplex[idx_b].i + cpx_matrix_a->arrayComplex[idx_a].i * cpx_matrix_b->arrayComplex[idx_b].r;
				}
				idx_c = row_i * cpx_matrix_c->column + column_j;
				cpx_matrix_c->arrayComplex[idx_c].r = temp_cpx.r;
				cpx_matrix_c->arrayComplex[idx_c].i = temp_cpx.i;
				temp_cpx.r = temp_cpx.i = 0;
			}
		}
	}
	return cpx_matrix_c;
}

matrix_state *matmul(matrix_state* matrix_a, matrix_state* matrix_b)
{
	matrix_state* matrix_c = (matrix_state*)malloc(sizeof(matrix_state));
	int row_i, column_j, ij;
	int idx_a, idx_b, idx_c;
	float temp;
	temp = 0.f;
	
	if (matrix_a->column != matrix_b->row)
	{
		puts("ERROE: An incompatable matrix!\n");
		return;
	}
	else
	{
		InitMatrix(matrix_c, matrix_a->row, matrix_b->column);
		for (row_i = 0; row_i < matrix_c->row; row_i++)
		{
			for (column_j = 0; column_j < matrix_c->column; column_j++)
			{
				for (ij = 0; ij < matrix_a->column; ij++)
				{
					idx_a = row_i * matrix_a->column + ij;
					idx_b = ij * matrix_b->column + column_j;
					temp += matrix_a->array[idx_a] * matrix_b->array[idx_b];
				}
				idx_c = row_i * matrix_c->column + column_j;
				matrix_c->array[idx_c] = temp;
				temp = 0;
			}
		}
	}
	return matrix_c;
}

complex_matrix_state* cpx_matrix_transpose(complex_matrix_state* cpx_matrix)
{
	complex_matrix_state* cpx_matrix_transpose = (complex_matrix_state*)malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(cpx_matrix_transpose, cpx_matrix->column, cpx_matrix->row);
	
	for (int i = 0; i < cpx_matrix->row; i++)
	{
		for (int j = 0; j < cpx_matrix->column; j++)
		{
			cpx_matrix_transpose->arrayComplex[j * cpx_matrix_transpose->column + i].r = cpx_matrix->arrayComplex[i * cpx_matrix->column + j].r;
			cpx_matrix_transpose->arrayComplex[j * cpx_matrix_transpose->column + i].i = -1 * cpx_matrix->arrayComplex[i * cpx_matrix->column + j].i;
		}
	}
	return cpx_matrix_transpose;
}


complex_matrix_state* cpx_matrix_inverse(complex_matrix_state* cpx_matrix)
{
	int row, column;
	int i, j, k, t;
	float norm;
	row = cpx_matrix->row;
	column = cpx_matrix->column;
	complex_matrix_state* L, *U, *L_Inverse, *U_Inverse;
	complex_matrix_state* cpx_matrix_inverse;
	kiss_fft_cpx temp;
	L = malloc(sizeof(complex_matrix_state));
	U = malloc(sizeof(complex_matrix_state));
	L_Inverse = malloc(sizeof(complex_matrix_state));
	U_Inverse = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(L, row, column);
	InitComplexMatrix(U, row, column);
	InitComplexMatrix(L_Inverse, row, column);
	InitComplexMatrix(U_Inverse, row, column);

	// 对角元素
	for (i = 0; i < row; i++)
	{
		L->arrayComplex[i * row + i].r = 1.f;
		L->arrayComplex[i * row + i].i = 0.f;
	}

	// U 的 第一行与matrix的第一行相等 
	for (j = 0; j < column; j++)
	{
		U->arrayComplex[j].r = cpx_matrix->arrayComplex[j].r;
		U->arrayComplex[j].i = cpx_matrix->arrayComplex[j].i;
	}

	// 比较系数
	norm = U->arrayComplex[0].r * U->arrayComplex[0].r + U->arrayComplex[0].i + U->arrayComplex[0].i;
	for (i = 1;i < column;i++)
	{
		L->arrayComplex[i * column].r = (cpx_matrix->arrayComplex[i * column].r * U->arrayComplex[0].r + cpx_matrix->arrayComplex[i * column].i * U->arrayComplex[0].i) / norm;
		L->arrayComplex[i * column].i = (cpx_matrix->arrayComplex[i * column].i * U->arrayComplex[0].r - cpx_matrix->arrayComplex[i * column].r * U->arrayComplex[0].i) / norm;
	}

	for (k = 1; k < row; k++)
	{
		for (j = k; j < column; j++)
		{
			temp.r = 0.0;
			temp.i = 0.0;
			for (t = 0; t < k; t++)
			{
				temp.r += L->arrayComplex[k * column + t].r * U->arrayComplex[t * column + j].r - L->arrayComplex[k * column + t].i * U->arrayComplex[t * column + j].i;
				temp.i += L->arrayComplex[k * column + t].r * U->arrayComplex[t * column + j].i + L->arrayComplex[k * column + t].i * U->arrayComplex[t * column + j].r;
			}
			U->arrayComplex[k * column + j].r = cpx_matrix->arrayComplex[k * column + j].r - temp.r;
			U->arrayComplex[k * column + j].i = cpx_matrix->arrayComplex[k * column + j].i - temp.i;
		}
		for (i = k; i < column; i++)
		{
			temp.r = 0.0;
			temp.i = 0.0;
			for (t = 0; t < k; t++)
			{
				temp.r += L->arrayComplex[i * column + t].r * U->arrayComplex[t * column + k].r - L->arrayComplex[i * column + t].i * U->arrayComplex[t * column + k].i;
				temp.i += L->arrayComplex[i * column + t].r * U->arrayComplex[t * column + k].i + L->arrayComplex[i * column + t].i * U->arrayComplex[t * column + k].r;
			}
			temp.r = cpx_matrix->arrayComplex[i * column + k].r - temp.r;
			temp.i = cpx_matrix->arrayComplex[i * column + k].i - temp.i;
			norm = U->arrayComplex[k * column + k].r * U->arrayComplex[k * column + k].r + U->arrayComplex[k * column + k].i * U->arrayComplex[k * column + k].i;
			L->arrayComplex[i * column + k].r = (temp.r * U->arrayComplex[i * column + k].r + U->arrayComplex[i * column + k].i * temp.i) / norm;
			L->arrayComplex[i * column + k].i = (temp.i * U->arrayComplex[i * column + k].r - U->arrayComplex[i * column + k].r * temp.i) / norm;
		}
	}

	for (i = 0; i < row; i++)
	{
		L_Inverse->arrayComplex[i * column + i].i = 0.0;
		L_Inverse->arrayComplex[i * column + i].r = 1.0;
	}

	for (j = 0; j < column; j++)
	{
		for (i = j + 1; i < row; i++)
		{
			temp.r = 0.0;
			temp.i = 0.0;
			for (k = j; k < i; k++) 
			{
				temp.r += L->arrayComplex[i * column + k].r * L_Inverse->arrayComplex[k * column + j].r - L->arrayComplex[i * column + k].i * L_Inverse->arrayComplex[k * column + j].i;
				temp.i += L->arrayComplex[i * column + k].i * L_Inverse->arrayComplex[k * column + j].r + L->arrayComplex[i * column + k].r * L_Inverse->arrayComplex[k * column + j].i;
			}
			L_Inverse->arrayComplex[i * column + j].r = (-1) * (L_Inverse->arrayComplex[j * column + j].r * temp.r - L_Inverse->arrayComplex[j * column + j].i * temp.i);
			L_Inverse->arrayComplex[i * column + j].i = (-1) * (L_Inverse->arrayComplex[j * column + j].r * temp.i + L_Inverse->arrayComplex[j * column + j].i * temp.r);
		}
	}
	kiss_fft_cpx temp2;
	temp2.r = 1.0;
	temp2.i = 0.0;
	for (i = 0; i < column; i++)                    
	{
		norm = U->arrayComplex[i * column + i].r * U->arrayComplex[i * column + i].r + U->arrayComplex[i * column + i].i * U->arrayComplex[i * column + i].i;
		U_Inverse->arrayComplex[i * column + i].r = (temp2.r * U->arrayComplex[i * column + i].r + temp2.i * U->arrayComplex[i * column + i].i) / norm;
		U_Inverse->arrayComplex[i * column + i].i = (temp2.i * U->arrayComplex[i * column + i].r - temp2.r * U->arrayComplex[i * column + i].i) / norm;
	}

	for (j = 0; j < column; j++)
	{
		for (i = j - 1; i >= 0; i--)
		{
			temp.r = 0.0;
			temp.i = 0.0;
			for (k = i + 1; k <= j; k++)
			{
				temp.r += U->arrayComplex[i * column + k].r * U_Inverse->arrayComplex[k * column + j].r - U->arrayComplex[i * column + k].i * U_Inverse->arrayComplex[k * column + j].i;
				temp.i += U->arrayComplex[i * column + k].i * U_Inverse->arrayComplex[k * column + j].r + U->arrayComplex[i * column + k].r * U_Inverse->arrayComplex[k * column + j].i;
			}
			temp.i = -temp.i;
			temp.r = -temp.r;

			norm = (U->arrayComplex[i * column + i].r * U->arrayComplex[i * column + i].r + U->arrayComplex[i * column + i].i + U->arrayComplex[i * column + i].i);
			U_Inverse->arrayComplex[i * column + j].r = (temp.r * U->arrayComplex[i * column + i].r + temp.i * U->arrayComplex[i * column + i].i) / norm;
			U_Inverse->arrayComplex[i * column + j].i = (temp.i * U->arrayComplex[i * column + i].r - temp.r * U->arrayComplex[i * column + i].i) / norm;
				
		}
	}
	cpx_matrix_inverse = cpx_matmul(U_Inverse, L_Inverse);
	return cpx_matrix_inverse;
}

void matrix_destroy(matrix_state* matrix)
{
	free(matrix->array);
	matrix->row = matrix->column = 0;
}

void cpx_matrix_destroy(complex_matrix_state* cpx_matrix)
{
	free(cpx_matrix->arrayComplex);
	cpx_matrix->row = cpx_matrix->column = 0;
}