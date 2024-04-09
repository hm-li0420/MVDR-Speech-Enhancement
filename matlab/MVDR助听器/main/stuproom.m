function y = stuproom(speech,cfg)
%生成多通道声音
%speech:原始语音
%c:声速 fs：采样频率 room：Room dimensions [x y z] (m) T60:时延

%% room setup
c = 343;                        % Sound velocity (m/s)
fs = cfg.fs;                    % Sample frequency (samples/s)
room = cfg.room;                % Room dimensions [x y z] (m)
T60 = cfg.T60;                          % Reverberation time (s)
beta = max(T60 * fs, 300);      % 混响点数
Nch = cfg.Nch;                              % Number of channels
micCenter = cfg.micCenter;          % array center (m)
micCoordinate= cfg.micCoordinate;   % microphone array coordinates
micPose = micCoordinate + repmat(micCenter, Nch, 1);

%% show microphone array
az = cfg.az;
el = cfg.el;
dist = cfg.dist;
sourcePose = dist * [cos(el/180*pi)*cos(az/180*pi) cos(el/180*pi)*sin(az/180*pi) sin(el/180*pi)] + micCenter;
%M * 3，M个mic，3维空间；

rirSimu = rir_generator(c, fs, micPose, sourcePose, room, T60, beta);%生成声音传播的空间滤波器（时域输出）

% observe
y = zeros(length(speech), Nch);
for ch=1:Nch
    y(:,ch) = fftfilt(rirSimu(ch,:), speech);
end
y = 0.8*y/max(max(abs(y)));
end