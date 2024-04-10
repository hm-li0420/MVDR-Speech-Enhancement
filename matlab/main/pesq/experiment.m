clear all; clc;
close all;
snr=[-5 0 5 10];
s1=[5/60*100 8/60*100 9/60*100;28/60*100 36/60*100 29/60*100;51/60*100 57/60*100 55/60*100;57/60*100 60/60*100 58/60*100];
figure('color','w');
b1=bar(s1,1);%1ÊÇwidth
set(gca, 'xticklabel', {'-5','0','5','10'});
legend('CIS','CIS-B','CIS-P');  
xlabel('Input  SNR(dB)');  
ylabel('Percentage   correct:%'); 
set(gca,'FontName','Times New Roman','FontSize',14);

s2=[4/60*100 5/60*100 6/60*100;25/60*100 41/60*100 35/60*100;54/60*100 56/60*100 55/60*100;58/60*100 60/60*100 60/60*100];
figure('color','w');
b2=bar(s2,1);%1ÊÇwidth
set(gca, 'xticklabel', {'-5','0','5','10'});im_hatchC = applyhatch_plusC(1,'+|-','krb');
% set(b2(1),'FaceColor','w');
% set(b2(2),'FaceColor',[1 0.5 1]);
% set(b2(3),'FaceColor',[1 1 0]);
% im_hatchC = applyhatch_plusC(1,'\-x.',[0 0 0;0.5 0 0;0.5 0.5 0.5;1 0 0]);
legend('CIS','CIS-B','CIS-P');  
xlabel('Input  SNR(dB)');  
ylabel('Percentage   correct:%'); 
set(gca,'FontName','Times New Roman','FontSize',14);

s3=[5/60*100 9/60*100 6/60*100;43/60*100 49/60*100 44/60*100;57/60*100 58/60*100 57/60*100;59/60*100 60/60*100 60/60*100];
figure('color','w');
b3=bar(s3,1);%1ÊÇwidth
set(gca, 'xticklabel', {'-5','0','5','10'});
legend('CIS','CIS-B','CIS-P');  
xlabel('Input  SNR(dB)');  
ylabel('Percentage   correct:%'); 
set(gca,'FontName','Times New Roman','FontSize',14);