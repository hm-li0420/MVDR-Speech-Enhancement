%% 初始化
clc;clear all;close all;
addpath('..\STFT\')                         %添加函数路径
addpath('..\Simulation\')
addpath('..\main\pesq')
addpath('..\Simulation\RIR-Generator\')
speechDir = '..\Simulation\Data\';          % 音频数据存储地址
saveDir = 'GeneratedDataULA\';              % 输出数据地址
if ~exist(saveDir)
    mkdir(saveDir)                          % 若无创建一个新的GeneratedDataULA文件夹
end
speechFile = 'fajw0_sa1.wav';               % 声源 speech.wav
intfFile = 'mrgg0_si1829.wav';              %干扰声
windNmae = 'wind_normal.wav';

[speech,fs] = audioread([speechDir speechFile]);    %读取音频
speech = speech(:,1);                               % 如果存储的是多声道音频，取第一声道
%% parameter 参数
M = 4    ;                               % 麦克风数量
d =  1.73  ;                                % 麦克风间距，单位cm
Nfft = 480;                               % STFT窗长度Nfft=128有更大提升
SNR = inf;                               %信噪比
SIR = inf;                                 %信干比  %若SNR / SIR=inf，则不加噪声 / 干扰
%% configuration 配置
c = 343;                                    %声速
cfg = [];
cfg.fs = fs;                         % sampling rate ，采样率一定要匹配
cfg.Nch = M;
cfg.c =c;
cfg.room = [12 6 6];                 % room dimension (m) 房间大小x,y.z
cfg.T60 = 0;                         % reverberation time (s) 混响RT60
cfg.micCenter = [6 3 3];             % array center (m)  %阵列中心
cfg.micCoordinate = [(0 : -0.01*d : -0.01*(M-1)*d)' zeros(cfg.Nch,2)]; %麦克风配置
cfg.Nfft = Nfft;
cfg.SNR = SNR;                % signal to noise rate
cfg.windFile = [speechDir windNmae];
cfg.SIR = SIR;                 % signal to interference rate
cfg.interfFile = [speechDir intfFile]; %干扰声路径参数

%声源空间位置
cfg.az =90;                         % azimuth水平角/度
cfg.el = 0;                         % elevation仰角、度
cfg.dist  = 2;                      %距离/m

cfg.azITF = 180;
cfg.elITF = 0;
cfg.distITF = 2;

%% 生成仿真的阵列信号
x = stuproom(speech,cfg);
y = add_intf_wind(x,cfg);
y = y-mean(y);

%%
Y = stft_multi_2(y, Nfft);
[fn,len,M]=size(Y);
%% 以上为生成仿真信号，可不看。

%% MVDR 参数
f = (0:Nfft/2)*fs/Nfft;
f(f==0)=0.001;      % w不能有0，否则后面1/w是无穷大
w = 2*pi*f;
tau0 = 0.01*d/c;        %τ0
Ytr = zeros(M,1);
Ytr2 = zeros(len,M);
w_mvdr = zeros(M,len);
R0 = zeros(M,M);
R = zeros(M,M);
RR = zeros(M,M);
%% DOA参数
theata = 0:5:180;           % 遍历角度
f0 = 1000;                  % 遍历频率
a=0;                      % 估计方向平滑参数

N =floor(f0*Nfft/fs)+1;     %计算f0所在频率点
f0 = (N-1)*fs/Nfft;

%% MVDR 
for ii = 1:fn
    %%find DOA
    Ytr(:,:) = Y(ii,N,:);           %len*M
    RR0 = Ytr*Ytr';                %M*M
    RR = 0.5*RR0+0.5*RR + eye(M);
%     if ii>1                         %平滑协方差矩阵
%         RR = 0.5*RR0+0.5*RR + eye(M);
%     else
%         RR=RR0+ eye(M) ;
%     end
    for t =1:length(theata)         %遍历
        dw = exp(-1j*2*pi*f0*0.01*d/c*cosd(theata(t))*[0:M-1]');
        P(ii,t) = 1/ abs(dw'* inv(RR) *dw);
    end
    [Pmax,index]=max( P(ii,:));         %谱搜索
    thetad(ii) = theata(index);
    if ii>2
        thetad(ii) = a*thetad(ii-1)+ (1-a)*thetad(ii);
    end
    %%    求解MVDR weight
    for n =1:round(Nfft/2+1)
        SV = exp(-1j*w(n)*[0:M-1]'*tau0*cosd(thetad(ii))) ;%thetad(ii)
        Ytr(:) = squeeze(Y(ii,n,:));
        R0 = Ytr * Ytr';
        ev = eig(R0);
%         ev = 1;
        if ii==1
            R = R0 +eye(M)*max(ev)/10;
        else
            R = ( R + R0 )/2 +eye(M)*max(ev)/100;   %(旧的协方差矩阵+新的协方差矩阵)/2    eye(M)对角加载，防止奇异矩阵。
        end
        invR = pinv(squeeze(R))*SV;
        t = invR /( SV'* invR );
        w_mvdr(:,n,ii) = invR /( SV'* invR );
    end
    %% processing start
    Ytr2(:,:) = squeeze(Y(ii,:,:));
    Xest(ii,:) = sum(  conj(w_mvdr(:,:,ii)) .* permute(Ytr2,[2,1]),1);  %滤波
    end
figure
plot(thetad)
xest = istft_multi_2(Xest, length(speech));


%% 保存音频
audiowrite([saveDir 'source.wav'],x(:,1),fs)
for i=1:M
    audiowrite([saveDir [num2str(i), 'noisy.wav']],y(1:length(speech),i),fs);
end
% audiowrite([saveDir 'noisy.wav'],y(1:length(speech),1),fs)
audiowrite([saveDir 'MVDR.wav' ] ,xest,fs);
%% omlsa
[in,out]=omlsa([saveDir 'noisy.wav'],[saveDir 'omlsa.wav']);
[in,out]=omlsa([saveDir 'MVDR.wav'],[saveDir 'MQomlsa.wav']);
%%
pesq0 = pesq([saveDir 'source.wav'],[saveDir 'noisy.wav'])

pesq_omlsa = pesq([saveDir 'source.wav'],[saveDir 'omlsa.wav'])
pesq_MVDR = pesq([saveDir 'source.wav'],[saveDir 'MVDR.wav'])
pesq_MQ = pesq([saveDir 'source.wav'],[saveDir 'MQomlsa.wav'])






