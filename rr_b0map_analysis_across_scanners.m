clear all, close all, clc
dirName = 'E:\scanner_data\twix_data\240731\';
fileName = 'meas_MID00079_FID281176_gre_b0map_1_5T_OFFRES';
% fileName = 'meas_MID00163_FID419272_gre_b0map_6_11_3T_OFFRES';
toImport = [dirName fileName '.mat'];
offres15 = load(toImport);

fileName = 'meas_MID00162_FID419271_gre_b0map_6_8_3T_OFFRES';
toImport = [dirName fileName '.mat'];
offres3 = load(toImport);

fileName = 'meas_MID00073_FID100464_b0_mapping_rudy_0_55T_OFFRES';
toImport = [dirName fileName '.mat'];
offres055 = load(toImport);

%%
figure, sliceViewer(offres055.LPF)
%%
idx = 28;

M = max(max([offres055.LPF(:,:,idx);offres15.LPF(:,:,idx);offres3.LPF(:,:,idx)]));
m = min(min([offres055.LPF(:,:,idx);offres15.LPF(:,:,idx);offres3.LPF(:,:,idx)]));

clims=[m M];
%%
figure, 
subplot(131), imagesc(offres055.LPF(:,:,idx),clims), colormap('gray'), colorbar
subplot(132), imagesc(offres15.LPF(:,:,idx),clims), colormap('gray'), colorbar
subplot(133), imagesc(offres3.LPF(:,:,idx),clims), colormap('gray'), colorbar


%%
imagine(cat(3,offres055.LPF(:,:,idx),offres15.LPF(:,:,idx),offres3.LPF(:,:,idx)))

%%
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(offres3.LPF(:,:,idx)));
axis image;
rois{1} = get_rois; %roi 0.55T
% %%
% figure;
% imagesc(abs(offres15.LPF(:,:,idx)));
% axis image;
% rois{2} = get_rois; %roi 1.5T
% %%
% figure;
% imagesc(abs(offres3.LPF(:,:,idx)));
% axis image;
% rois{3} = get_rois; %roi 3T
%%

rois{2} = rois{1};
rois{3} = rois{1};

%%
[ix,iy] = find(rois{1}{1}==1);
fp055 = reshape(offres055.LPF(ix,iy,idx),[size(ix,1)*size(iy,1),1]);
avg055 = mean(fp055);
std055 = std(fp055);
disp(['0.55T: ' num2str(avg055) ' /pm ' num2str(std055) ' Hz' ])

[ix,iy] = find(rois{2}{1}==1);
fp15 = reshape(offres15.LPF(ix,iy,idx),[size(ix,1)*size(iy,1),1]);
avg15 = mean(fp15);
std15 = std(fp15);
disp(['1.5T: ' num2str(avg15) ' /pm ' num2str(std15) ' Hz' ])

[ix,iy] = find(rois{3}{1}==1);
fp3 = reshape(offres3.LPF(ix,iy,idx),[size(ix,1)*size(iy,1),1]);
avg3 = mean(fp3);
std3 = std(fp3);
disp(['3T: ' num2str(avg3) ' /pm ' num2str(std3) ' Hz' ])


