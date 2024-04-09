#ifndef MVDR_H
#define MVDR_H
#include "kiss_fftr.h"
#include "matrix.h"

#define FRAMESIZE 480
#define STEPSIZE 240
#define MIC_NUMS 4
#define MIC_D 0.0173
#define SAMPLERATES 16000

#define SEARCH_ANGLE_LOW 0
#define SEARCH_ANGLE_HIGH 180
#define SEARCH_ANGLE_STEP 5
#define SEARCH_FREQ 1000
#define ALPHA 0.5


typedef struct BF_State 
{
	kiss_fft_scalar* mic0;
	kiss_fft_scalar* mic1;
	kiss_fft_scalar* mic2;
	kiss_fft_scalar* mic3;

	kiss_fft_cpx* mic0_cout;
	kiss_fft_cpx* mic1_cout;
	kiss_fft_cpx* mic2_cout;
	kiss_fft_cpx* mic3_cout;
	kiss_fft_cpx* output;

	complex_matrix_state* guide_vector;
	complex_matrix_state* guide_vector_T;

	complex_matrix_state* Ryy_inverse;
	complex_matrix_state* Ryy;
	complex_matrix_state* R;

	complex_matrix_state* Ytr;
	complex_matrix_state* Ytr_T;

	complex_matrix_state* Dw;
	complex_matrix_state* Dw_T;

	complex_matrix_state* tempa;
	complex_matrix_state* tempb;
	complex_matrix_state* mvdr_a;
	complex_matrix_state* mvdr_c;
	kiss_fft_cpx* mvdr_weights;

	kiss_fftr_cfg kiss_fftr_state;
	kiss_fftr_cfg kiss_ifftr_state;
}bf_state;

bf_state* mvdr_init(int frame_size, int mic_nums);
void beamform(short* pcm_buff, int pcm_len, bf_state* bf_st, float angle);
int get_angle(short* pcm_buff, int pcm_len, bf_state* bf_st, int low_angle, int high_angle, int step_angle, float search_freq, float alpha);
void bf_destory(bf_state* st);
#endif