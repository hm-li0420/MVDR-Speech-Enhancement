function y = add_intf_wind(source,cfg)
%%
%  加干扰声和噪声
%%
SIR = cfg.SIR;
c = cfg.c; 
fs = cfg.fs;
room = cfg.room; 
T60 = cfg.T60;
Nch =cfg.Nch;
micCenter = cfg.micCenter;
micCoordinate= cfg.micCoordinate;
micPose = micCoordinate + repmat(micCenter, Nch, 1);
beta = max(T60 * fs, 300);      % 混响点数
%% add intf
y = source;                     %初始化
if SIR ~= inf
    %干扰信号空间位置信息
    azITF = cfg.azITF; elITF = cfg.elITF; distITF = cfg.distITF;
    interfPose = micCenter + distITF * [cos(elITF/180*pi)*cos(azITF/180*pi) cos(elITF/180*pi)*sin(azITF/180*pi) sin(elITF/180*pi)];
   %生成干扰信号
    rirITFSimu = rir_generator(c, fs, micPose, interfPose, room, T60, beta);
    %读取干扰信号
    interfFile = cfg.interfFile;
    interfSource = audioread([interfFile]);
    interfSource = interfSource(:,1);
    if size(interfSource,1) > size(source,1)  %截取与原始语音相同长度干扰语音
%         fst = floor(rand * (size(interfSource,1) - size(source,1)));   %随机取两段话长度差的起点
        fst=0;
        interfSource = interfSource(fst+1:fst+size(source,1));
    elseif size(interfSource,1) < size(source,1)            %过短补零
        interfSource = [interfSource;zeros(size(source,1)-size(interfSource,1),1)];
    end 
    interf = zeros(length(source), Nch);
    for ch=1:Nch
        interf(:,ch) =  fftfilt(rirITFSimu(ch,:), interfSource);
    end
    
    for n =1:Nch  %每个通道加 SIR干扰音    
        interf(:,n) = sqrt( sum(source(:,n).^2) / sum(interf(:,n).^2)/(10^(SIR/10)) ) * interf(:,n);
        y(:,n) = y(:,n) + interf(:,n);
    end
end
%% add noise
SNR =cfg.SNR ;

if SNR ~= inf
    [wind,fs] = audioread(cfg.windFile);
    wind = wind(:,1);
    for n =1:Nch  %每个通道加 SNR干扰音
        if size(wind,1) >= size(source,1)  %截取与原始语音相同长度干扰语音
            fst = floor(rand * (size(wind,1) - size(source,1)));   %随机取两段话长度差的起点
            wind_c = wind(fst+1:fst+size(source,1));
        else             %过短补零
            wind_c = [wind;zeros(size(source,1)-size(wind,1),1)];
        end
        wind_c(:) = sqrt( sum(source(:,n).^2) / sum(wind_c(:).^2)/(10^(SNR/10)) ) * wind_c(:);
        
        y(:,n) = y(:,n) + wind_c(:);      %加入风噪
    end
end
end