clear all, close all, clc
folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
fatxtfile = [folderTXT '\FISP_FA_body.txt'];

log.FA = importdata(fatxtfile);
Nex = length(log.FA); %1000;
log.FA15 = round(log.FA*1.5,2);
log.FA2 = round(log.FA*2,2);

% plot the parameters
figure('name','Acquisition Parameters');
plot(log.FA(1:Nex,1),'LineWidth',2); hold on
plot(log.FA15(1:Nex,1),'LineWidth',2); hold on
plot(log.FA2(1:Nex,1),'LineWidth',2); hold on
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');

%% translation to have Nex 400 or 300 
s1 = log.FA(210:end);
s2 = log.FA(1:209);
log.FA_shift = [s1;s2];

%%
figure ,plot(log.FA), hold on
plot(log.FA_shift,'r')
%%
savedir = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp\';
fileID = fopen([savedir 'FISP_FA_body_shifted.txt'],'w');
formatSpec = '%.2f\n';
for row = 1:Nex
    fprintf(fileID,formatSpec,log.FA_shift(row));
end
fclose(fileID);
