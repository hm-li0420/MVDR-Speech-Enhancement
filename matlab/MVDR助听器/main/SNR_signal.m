function snr=SNR_signal(I,In)
% ������������źŵ������
% I �Ǵ������ź�
% In �Ǵ���������ź�
% ����ȼ��㹫ʽ��
% snr=10*log10(Esignal/Enoise)

Ps=sum((I).^2);              % �źŵ�����
Pn=sum((I-In).^2);                   % ����������
snr=10*log10(Ps/Pn);                 % �źŵ�����������������֮�ȣ�����ֱ�ֵ