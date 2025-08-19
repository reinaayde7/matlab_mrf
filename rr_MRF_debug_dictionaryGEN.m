clear all, close all, clc
%%
isTRUEFISP = 1;
disp('open Siemens Raw data')
RawData.folder = 'E:\scanner_data\twix_data\240628\';

if isTRUEFISP
    RawData.name = 'meas_MID00286_FID96913_rrMRFv2_trueFISP.dat';
else
    RawData.name = 'meas_MID00285_FID96912_rrMRFv2_FISP.dat';
end
RawDataFileName = fullfile(RawData.folder,RawData.name);

[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,rawinfo,~] = loadSiemensRawData(RawDataFileName);
size(raw)
%%

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
flipangle = importdata([folderTXT '\FISP_FA.txt']);
tr0 = rawinfo.TR(1)/1000/1000;
tr = importdata([folderTXT '\FISP_TR.txt'])*1e-6 + tr0;

% plot the parameters
figure('name','Acquisition Parameters');
subplot(211);plot(flipangle(1:rawinfo.Nex,1),'LineWidth',2);
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');
subplot(212);plot(tr(1:rawinfo.Nex,1)*1000,'LineWidth',2);
title('Repetition Time'); xlabel('Time Points');ylabel('TR (ms)');


%% Calculate dictionary based on input text file

delay = 3000;
t1series = 1000; %[10:10:2000 2020:20:3000 3050:50:3500 4000:500:5000];
t2series = 200; %[2:2:100 105:5:300 310:10:500 520:20:800];% 850:50:1200];% 310:10:500 520:20:800 850:50:1500 1600:100:2000];

trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
fatxtfile = [folderTXT '\FISP_FA.txt'];
phtxtfile = [folderTXT '\FISP_PH.txt'];


if isTRUEFISP
    [dict,r] = Calculate_MRF_TRUEFISP_DictwithDelays(RawDataFileName,rawinfo,t1series,t2series,0,...%[0.6:0.1:1.4] for B1
    tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,delay,0);
    % save('Dict_TRUEFISP_noSPcorr_v5.mat','dict','r','-v7.3');
else
    [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,rawinfo,t1series,t2series,0,...%[0.6:0.1:1.4] for B1
    tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,delay,2);
    % save('Dict_FISP_noSPcorr.mat','dict','r','-v7.3');
end
  

figure,
subplot(311), plot(abs(dict)), ylabel('abs')
subplot(312), plot(real(dict)), ylabel('real')
subplot(313), plot(imag(dict)), ylabel('imag')



%%
figure,
subplot(311), plot(abs(squeeze(image_combined(200,200,:)))), ylabel('abs')
subplot(312), plot(real(squeeze(image_combined(200,200,:)))), ylabel('real')
subplot(313), plot(imag(squeeze(image_combined(200,200,:)))), ylabel('imag')