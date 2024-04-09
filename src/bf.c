#include "bf.h"
#include "kiss_fftr.h"

#define M_PI 3.14159265358979323846
extern short* mic0_input_buffer, * mic1_input_buffer, * mic2_input_buffer, * mic3_input_buffer;
extern float* mic0_out, * mic1_out, * mic2_out, * mic3_out;
extern float* mic0_pre;
float* windows;

bf_state* mvdr_init(int frame_size, int mic_nums)
{	
	mic0_input_buffer = malloc(frame_size * sizeof(short));
	mic1_input_buffer = malloc(frame_size * sizeof(short));
	mic2_input_buffer = malloc(frame_size * sizeof(short));
	mic3_input_buffer = malloc(frame_size * sizeof(short));

	mic0_out = malloc(frame_size * sizeof(float));
	mic1_out = malloc(frame_size * sizeof(float));
	mic2_out = malloc(frame_size * sizeof(float));
	mic3_out = malloc(frame_size * sizeof(float));

	mic0_pre = malloc(frame_size / 2 * sizeof(float));

	memset(mic0_input_buffer, 0, frame_size * sizeof(short));
	memset(mic1_input_buffer, 0, frame_size * sizeof(short));
	memset(mic2_input_buffer, 0, frame_size * sizeof(short));
	memset(mic3_input_buffer, 0, frame_size * sizeof(short));

	memset(mic0_pre, 0, frame_size / 2 * sizeof(float));

	bf_state* st = (bf_state*)malloc(sizeof(bf_state));

	st->mic0 = malloc(frame_size * sizeof(kiss_fft_scalar));
	st->mic1 = malloc(frame_size * sizeof(kiss_fft_scalar));
	st->mic2 = malloc(frame_size * sizeof(kiss_fft_scalar));
	st->mic3 = malloc(frame_size * sizeof(kiss_fft_scalar));

	st->mic0_cout = malloc(frame_size * sizeof(kiss_fft_cpx));
	st->mic1_cout = malloc(frame_size * sizeof(kiss_fft_cpx));
	st->mic2_cout = malloc(frame_size * sizeof(kiss_fft_cpx));
	st->mic3_cout = malloc(frame_size * sizeof(kiss_fft_cpx));
	st->output = malloc(frame_size * sizeof(kiss_fft_cpx));

	st->guide_vector = (complex_matrix_state*)malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->guide_vector, mic_nums, 1);

	st->guide_vector_T = (complex_matrix_state*)malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->guide_vector_T, 1, mic_nums);

	st->Ryy_inverse = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Ryy_inverse, mic_nums, mic_nums);

	st->Ryy = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Ryy, mic_nums, mic_nums);

	st->R = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->R, mic_nums, mic_nums);

	st->Ytr = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Ytr, mic_nums, 1);

	st->Ytr_T = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Ytr_T, 1, mic_nums);

	st->Dw = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Dw, mic_nums, 1);

	st->Dw_T = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->Dw_T, 1, mic_nums);

	st->tempa = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->tempa, 1, mic_nums);

	st->tempb = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->tempb, 1, 1);

	st->mvdr_a = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->mvdr_a,  mic_nums, 1);

	st->mvdr_c = malloc(sizeof(complex_matrix_state));
	InitComplexMatrix(st->mvdr_c, 1, 1);

	st->mvdr_weights = malloc(mic_nums * sizeof(kiss_fft_cpx));

	st->kiss_fftr_state = kiss_fftr_alloc(frame_size, 0, 0, 0);
	st->kiss_ifftr_state = kiss_fftr_alloc(frame_size, 1, 0, 0);


	memset(st->mic0, 0, frame_size * sizeof(kiss_fft_scalar));
	memset(st->mic1, 0, frame_size * sizeof(kiss_fft_scalar));
	memset(st->mic2, 0, frame_size * sizeof(kiss_fft_scalar));
	memset(st->mic3, 0, frame_size * sizeof(kiss_fft_scalar));
	memset(st->output, 0, frame_size * sizeof(kiss_fft_scalar));

	windows = malloc(frame_size * sizeof(float));
	for (int i = 0; i < frame_size; i++) {
		windows[i] = 0.5f - 0.5f * cosf((float)(M_PI * 2 * ((float)i) / frame_size));
	}

	return st;
}

void beamform(short* pcm_buff, int pcm_len, bf_state* bf_st, float angle)
{
	float norm;
	int frame_size = 2 * pcm_len / MIC_NUMS;
	float tau = MIC_D / 340;
	float step_freq = SAMPLERATES / FRAMESIZE; 
	
	for (int i = 0; i < frame_size; i++)
	{
		if (i < frame_size / 2 + 1)
		{
			for (int j = 0;j < MIC_NUMS;j++)
			{
				bf_st->guide_vector->arrayComplex[j].r = cos(-1 * 2 * M_PI * i * step_freq * tau * j * cos(angle / 180.f * M_PI));
				bf_st->guide_vector->arrayComplex[j].i = sin(-1 * 2 * M_PI * i * step_freq * tau * j * cos(angle / 180.f * M_PI));
			}
			bf_st->guide_vector_T = cpx_matrix_transpose(bf_st->guide_vector);

			// 4 mics
			bf_st->Ytr->arrayComplex[0].r = bf_st->mic0_cout[i].r;
			bf_st->Ytr->arrayComplex[0].i = bf_st->mic0_cout[i].i;

			bf_st->Ytr->arrayComplex[1].r = bf_st->mic1_cout[i].r;
			bf_st->Ytr->arrayComplex[1].i = bf_st->mic1_cout[i].i;

			bf_st->Ytr->arrayComplex[2].r = bf_st->mic2_cout[i].r;
			bf_st->Ytr->arrayComplex[2].i = bf_st->mic2_cout[i].i;

			bf_st->Ytr->arrayComplex[3].r = bf_st->mic3_cout[i].r;
			bf_st->Ytr->arrayComplex[3].i = bf_st->mic3_cout[i].i;

			bf_st->Ytr_T = cpx_matrix_transpose(bf_st->Ytr);
			bf_st->R = cpx_matmul(bf_st->Ytr, bf_st->Ytr_T);

			for (int k = 0;k < MIC_NUMS;k++)
			{
				for (int q = 0;q < MIC_NUMS;q++)
				{
					if (k == q)
					{
						bf_st->R->arrayComplex[k * MIC_NUMS + q].r = bf_st->R->arrayComplex[k * MIC_NUMS + q].r + 1;
					}
				}
			}

			bf_st->R = cpx_matrix_inverse(bf_st->R);

			bf_st->mvdr_a = cpx_matmul(bf_st->R, bf_st->guide_vector);
			bf_st->mvdr_c = cpx_matmul(bf_st->guide_vector_T, bf_st->mvdr_a);

			for (int k = 0;k < MIC_NUMS;k++)
			{
				norm = bf_st->mvdr_c->arrayComplex->r * bf_st->mvdr_c->arrayComplex->r + bf_st->mvdr_c->arrayComplex->i * bf_st->mvdr_c->arrayComplex->i;
				bf_st->mvdr_weights[k].r = (bf_st->mvdr_a->arrayComplex[k].r * bf_st->mvdr_c->arrayComplex[0].r + bf_st->mvdr_a->arrayComplex[k].i * bf_st->mvdr_c->arrayComplex[0].i) / (norm + 1e-9);
				bf_st->mvdr_weights[k].i = (bf_st->mvdr_a->arrayComplex[k].i * bf_st->mvdr_c->arrayComplex[0].r - bf_st->mvdr_a->arrayComplex[k].r * bf_st->mvdr_c->arrayComplex[0].i) / (norm + 1e-9);

				bf_st->output[i].r += (bf_st->mvdr_weights[k].r * bf_st->Ytr->arrayComplex[k].r - bf_st->mvdr_weights[k].i * bf_st->Ytr->arrayComplex[k].i) / frame_size;
				bf_st->output[i].i += (bf_st->mvdr_weights[k].i * bf_st->Ytr->arrayComplex[k].r + bf_st->mvdr_weights[k].r * bf_st->Ytr->arrayComplex[k].i) / frame_size;
			}
		}
		else
		{
			bf_st->output[i].r = bf_st->output[frame_size - i].r;
			bf_st->output[i].i = -1*bf_st->output[frame_size - i].i;
		}
	}

	kiss_fftri(bf_st->kiss_ifftr_state, bf_st->output, bf_st->mic0);
	memset(bf_st->output, 0, sizeof(kiss_fft_cpx) * frame_size);
	
	for (int k = 0; k < frame_size; k++)
	{
		mic0_out[k] = bf_st->mic0[k];
	}

	// 重叠相加还原为一帧语音
	for (int i = 0; i < STEPSIZE; i++)
	{
		if ((mic0_pre[i] + mic0_out[i]) > 0)
		{
			pcm_buff[i] = (short)((mic0_pre[i] + mic0_out[i]) * 32767);
		}
		else
		{
			pcm_buff[i] = (short)((mic0_pre[i] + mic0_out[i]) * 32768);
		}
		if (pcm_buff[i] > 32767)
		{
			pcm_buff[i] = 32767;
		}
		else if (pcm_buff[i] < -32768)
		{
			pcm_buff[i] = -32768;
		}
	}

	for (int i = 0;i < frame_size - STEPSIZE; i++)
	{
		mic0_pre[i] = mic0_out[i + STEPSIZE];
	}
}


int get_angle(short* pcm_buff, int pcm_len, bf_state* bf_st, int low_angle, int high_angle, int step_angle, float search_freq, float alpha)
{
	float abs_p, abs_max;
	int angle = low_angle;
	int angle_max;
	int N = (int)((search_freq * FRAMESIZE / SAMPLERATES));
	int frame_size = 2 * pcm_len / MIC_NUMS;
	int step = frame_size / 2;
	float tau = MIC_D / 340;
	float step_freq = SAMPLERATES/FRAMESIZE;

	memcpy(mic0_input_buffer, &mic0_input_buffer[step], sizeof(short) * (frame_size - step));
	memcpy(mic1_input_buffer, &mic1_input_buffer[step], sizeof(short) * (frame_size - step));
	memcpy(mic2_input_buffer, &mic2_input_buffer[step], sizeof(short) * (frame_size - step));
	memcpy(mic3_input_buffer, &mic3_input_buffer[step], sizeof(short) * (frame_size - step));

	for (int i = 0;i < pcm_len / MIC_NUMS;i++)
	{
		mic0_input_buffer[step + i] = pcm_buff[i * MIC_NUMS];
		mic1_input_buffer[step + i] = pcm_buff[i * MIC_NUMS + 1];
		mic2_input_buffer[step + i] = pcm_buff[i * MIC_NUMS + 2];
		mic3_input_buffer[step + i] = pcm_buff[i * MIC_NUMS + 3];
	}

	// 语音加窗
	for (int i = 0; i < frame_size; i++)
	{
		if (mic0_input_buffer[i] > 0)
		{
			mic0_out[i] = (float)mic0_input_buffer[i] / 32767 * windows[i];
		}
		else
		{
			mic0_out[i] = (float)mic0_input_buffer[i] / 32768 * windows[i];
		}

		if (mic1_input_buffer[i] > 0)
		{
			mic1_out[i] = (float)mic1_input_buffer[i] / 32767 * windows[i];
		}
		else
		{
			mic1_out[i] = (float)mic1_input_buffer[i] / 32768 * windows[i];
		}

		if (mic2_input_buffer[i] > 0)
		{
			mic2_out[i] = (float)mic2_input_buffer[i] / 32767 * windows[i];
		}
		else
		{
			mic2_out[i] = (float)mic2_input_buffer[i] / 32768 * windows[i];
		}

		if (mic3_input_buffer[i] > 0)
		{
			mic3_out[i] = (float)mic3_input_buffer[i] / 32767 * windows[i];
		}
		else
		{
			mic3_out[i] = (float)mic3_input_buffer[i] / 32768 * windows[i];
		}
	}

	// FFT
	for (int i = 0; i < frame_size; ++i)
	{
		bf_st->mic0[i] = mic0_out[i];
		bf_st->mic1[i] = mic1_out[i];
		bf_st->mic2[i] = mic2_out[i];
		bf_st->mic3[i] = mic3_out[i];
	}

	kiss_fftr(bf_st->kiss_fftr_state, bf_st->mic0, bf_st->mic0_cout);
	kiss_fftr(bf_st->kiss_fftr_state, bf_st->mic1, bf_st->mic1_cout);
	kiss_fftr(bf_st->kiss_fftr_state, bf_st->mic2, bf_st->mic2_cout);
	kiss_fftr(bf_st->kiss_fftr_state, bf_st->mic3, bf_st->mic3_cout);

	bf_st->Ytr->arrayComplex[0].r = bf_st->mic0_cout[N].r;
	bf_st->Ytr->arrayComplex[0].i = bf_st->mic0_cout[N].i;

	bf_st->Ytr->arrayComplex[1].r = bf_st->mic1_cout[N].r;
	bf_st->Ytr->arrayComplex[1].i = bf_st->mic1_cout[N].i;

	bf_st->Ytr->arrayComplex[2].r = bf_st->mic2_cout[N].r;
	bf_st->Ytr->arrayComplex[2].i = bf_st->mic2_cout[N].i;

	bf_st->Ytr->arrayComplex[3].r = bf_st->mic3_cout[N].r;
	bf_st->Ytr->arrayComplex[3].i = bf_st->mic3_cout[N].i;

	bf_st->Ytr_T = cpx_matrix_transpose(bf_st->Ytr);
	bf_st->R = cpx_matmul(bf_st->Ytr, bf_st->Ytr_T);


	for (int i = 0;i < MIC_NUMS;i++) 
	{
		for (int j = 0;j < MIC_NUMS;j++)
		{
			if (i == j)
			{
				bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].r = ALPHA * bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].r + (1 - ALPHA) * bf_st->R->arrayComplex[i * MIC_NUMS + j].r + 1;
				bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].i = ALPHA * bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].i + (1 - ALPHA) * bf_st->R->arrayComplex[i * MIC_NUMS + j].i;
			}
			else
			{
				bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].r = ALPHA * bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].r + (1 - ALPHA) * bf_st->R->arrayComplex[i * MIC_NUMS + j].r;
				bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].i = ALPHA * bf_st->Ryy->arrayComplex[i * MIC_NUMS + j].i + (1 - ALPHA) * bf_st->R->arrayComplex[i * MIC_NUMS + j].i;
			}
		}
	}

	bf_st->Ryy_inverse = cpx_matrix_inverse(bf_st->Ryy);

	abs_max = 0.f;
	angle_max = 0;
	while(angle <= high_angle)
	{
		for (int i = 0;i < MIC_NUMS;i++)
		{
			bf_st->Dw->arrayComplex[i].r = cos(-1 * 2 * M_PI * search_freq * tau * cosf(angle / 180.f * M_PI) * i);
			bf_st->Dw->arrayComplex[i].i = sin(-1 * 2 * M_PI * search_freq * tau * cosf(angle / 180.f * M_PI) * i);
		}
		bf_st->Dw_T = cpx_matrix_transpose(bf_st->Dw);
		bf_st->tempa = cpx_matmul(bf_st->Dw_T, bf_st->Ryy_inverse);
		bf_st->tempb = cpx_matmul(bf_st->tempa, bf_st->Dw);
		abs_p = 1 / (bf_st->tempb->arrayComplex->r * bf_st->tempb->arrayComplex->r + bf_st->tempb->arrayComplex->i * bf_st->tempb->arrayComplex->i);
		if (abs_p > abs_max)
		{
			abs_max = abs_p;
			angle_max = angle;
		}
		angle += step_angle;
	}
	//printf("%d \n", angle_max);
	return angle_max;
}

void bf_destory(bf_state* st)
{
	free(st->mic0);
	free(st->mic1);
	free(st->mic2);
	free(st->mic3);

	free(st->mic0_cout);
	free(st->mic1_cout);
	free(st->mic2_cout);
	free(st->mic3_cout);

	free(st->output);
	
	cpx_matrix_destroy(st->guide_vector);
	cpx_matrix_destroy(st->guide_vector_T);

	cpx_matrix_destroy(st->Ryy_inverse);
	cpx_matrix_destroy(st->Ryy);
	cpx_matrix_destroy(st->R);
	
	cpx_matrix_destroy(st->Ytr);
	cpx_matrix_destroy(st->Ytr_T);

	cpx_matrix_destroy(st->Dw);
	cpx_matrix_destroy(st->Dw_T);

	cpx_matrix_destroy(st->tempa);
	cpx_matrix_destroy(st->tempb);
	cpx_matrix_destroy(st->mvdr_a);
	cpx_matrix_destroy(st->mvdr_c);
	free(st->mvdr_weights);

	kiss_fftr_free(st->kiss_fftr_state);
	kiss_fftr_free(st->kiss_ifftr_state);

	free(st);

	free(mic0_input_buffer);
	free(mic1_input_buffer);
	free(mic2_input_buffer);
	free(mic3_input_buffer);

	free(mic0_out);
	free(mic1_out);
	free(mic2_out);
	free(mic3_out);

	free(mic0_pre);
	free(windows);
}