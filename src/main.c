#include <stdio.h>
#include <stdint.h>
#include "matrix.h"
#include "bf.h"


short *mic0_input_buffer , *mic1_input_buffer, * mic2_input_buffer, * mic3_input_buffer;
float *mic0_out, *mic1_out, * mic2_out, * mic3_out;
float *mic0_pre;

 
int main(int argc, char **argv)
{
	int angle;
	if (argc != 6)
	{
		fprintf(stderr, "BF mic_0.pcm mic_1.pcm mic2.pcm mic3.pcm out.pcm\n");
		exit(1);
	}
	FILE* mic0_input_file = fopen(argv[1], "rb");
	if (mic0_input_file == NULL) {
		printf("Cannot open mic0.pcm \n");
		return -1;
	}

	FILE* mic1_input_file = fopen(argv[2], "rb");
	if (mic0_input_file == NULL) {
		printf("Cannot open mic1.pcm \n");
		return -1;
	}

	FILE* mic2_input_file = fopen(argv[3], "rb");
	if (mic0_input_file == NULL) {
		printf("Cannot open mic2.pcm \n");
		return -1;
	}

	FILE* mic3_input_file = fopen(argv[4], "rb");
	if (mic0_input_file == NULL) {
		printf("Cannot open mic3.pcm \n");
		return -1;
	}

	FILE* mic0_output_file = fopen(argv[5], "wb");
	if (mic0_input_file == NULL) {
		printf("Cannot open out.pcm \n");
		return -1;
	}

	bf_state* bf_st = mvdr_init(FRAMESIZE, MIC_NUMS);

	//定义PC仿真用的外部buffer
	int16_t mic0_external_buff[STEPSIZE];
	int16_t mic1_external_buff[STEPSIZE];
	int16_t mic2_external_buff[STEPSIZE];
	int16_t mic3_external_buff[STEPSIZE];
	int16_t pcm_external_buff[STEPSIZE * MIC_NUMS];
	int _pcm_len = STEPSIZE * MIC_NUMS;

	while (!feof(mic0_input_file) && !feof(mic1_input_file) && !feof(mic2_input_file) && !feof(mic3_input_file))
	{
		int pcm_len = _pcm_len;
		fread(mic0_external_buff, sizeof(short), STEPSIZE, mic0_input_file);
		fread(mic1_external_buff, sizeof(short), STEPSIZE, mic1_input_file);
		fread(mic2_external_buff, sizeof(short), STEPSIZE, mic2_input_file);
		fread(mic3_external_buff, sizeof(short), STEPSIZE, mic3_input_file);

		//将4 mic的输入信号依次间隔放入pcm_buff, 模拟BES平台的处理输入信号
		for (int i = 0; i < STEPSIZE; i++) 
		{
			pcm_external_buff[i * MIC_NUMS] = mic0_external_buff[i];
			pcm_external_buff[i * MIC_NUMS + 1] = mic1_external_buff[i];
			pcm_external_buff[i * MIC_NUMS + 2] = mic2_external_buff[i];
			pcm_external_buff[i * MIC_NUMS + 3] = mic3_external_buff[i];
		}

		// 麦克风阵列函数，输入为两个通道的语音,输出为单通道，_pcm_len 应为每个通道样点数*通道数,参考BES平台的定义
		angle = get_angle(pcm_external_buff, pcm_len, bf_st, SEARCH_ANGLE_LOW, SEARCH_ANGLE_HIGH, SEARCH_ANGLE_STEP, SEARCH_FREQ, ALPHA);
		beamform(pcm_external_buff, pcm_len, bf_st, angle);
		pcm_len /= MIC_NUMS;

		// 数据写入PCM文件
		fwrite(pcm_external_buff, sizeof(short), STEPSIZE, mic0_output_file);
	}
	bf_destory(bf_st);
	fclose(mic0_input_file);
	fclose(mic1_input_file);
	fclose(mic2_input_file);
	fclose(mic3_input_file);
	fclose(mic0_output_file);
	return 0;
}

