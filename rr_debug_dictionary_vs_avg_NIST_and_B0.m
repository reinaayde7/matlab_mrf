clear all, close all, clc

%% load measurments

basedir = 'E:\scanner_data\twix_data\241011\';

filename = 'avgMEAS_17-Oct-2024_meas_MID01026_fp_b0_nist';
load([basedir filename]); %variables: nist, meas


%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
% fatxtfile = [folderTXT '\FISP_FA_orig.txt'];
fatxtfile = [folderTXT '\FISP_FA_body.txt'];
phtxtfile = [folderTXT '\FISP_PH.txt'];

TR0 = 18000;
TI = 20640;
Nex = 1000;



log.FA = importdata(fatxtfile);
log.tr0 = TR0/1000/1000;
log.ti = TI/1000/1000;
log.TR = importdata(trtxtfile)*1e-6 + log.tr0;
log.TE = importdata(tetxtfile);
log.PH = importdata(phtxtfile);

log.FA=log.FA(1:Nex);
log.TR=log.FA(1:Nex);
log.TE=log.FA(1:Nex);
log.PH=log.FA(1:Nex);

% plot the parameters
figure('name','Acquisition Parameters');
subplot(211);plot(log.FA,'LineWidth',2);
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');
subplot(212);plot(log.TR*1000,'LineWidth',2);
title('Repetition Time'); xlabel('Time Points');ylabel('TR (ms)');

%%
vialN = 8;
log.delay = 3000;
log.t1series = nist.values(vialN,1); 
log.t2series = nist.values(vialN,2); 
log.offres = meas.avgB0(vialN); 

[dict,r] = DebugCalculate_MRF_TRUEFISP_DictwithDelays(log,log.t1series,log.t2series,log.offres,...%[0.6:0.1:1.4] for B1
tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);

save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_trueFISP_t1t2B0v13.mat'],'dict','r','log','-v7.3');

%%
