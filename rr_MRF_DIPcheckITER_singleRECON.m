% Rudy December 2024
% check DIP recon for SINGLE recon slice across saved iteration for quality on maps
clear all, close all, clc

%%
MRFdir = ['E:\POSTDOC_UoM\08_Project_MRF\DIP\TRUEFISP\Brain_HV4_250228_tilt_f300_Nex600_FAbody15\'];
MRFDIPname = ['Brain_HV4_250228_tilt_s20_f300_Nex600_FAbody15'];

%%
dataDIP = load([MRFdir MRFDIPname '\bSSFP_DIP_Drop5_Epochs400.mat']);

t1maps = rot90(squeeze(permute(dataDIP.T1Iter,[2,3,1])));
t2maps = rot90(squeeze(permute(dataDIP.T2Iter,[2,3,1])));
m0maps = rot90(squeeze(permute(dataDIP.M0Iter,[2,3,1])));

%%
b0maps = rot90(squeeze(permute(dataDIP.B0Iter,[2,3,1])));

%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 250]; %uplimit:1200 for NIST T2 at 0.55T 

toSliceViewer = t1maps;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = t2maps; 
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar
% 
toSliceViewer = abs(m0maps);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

%% 
toSliceViewer = b0maps;
figure,
sliceViewer(toSliceViewer, "DisplayRange", [-100, +100], 'ScaleFactors', [1,1,1]), colorbar

%%
figure, plot(dataDIP.lossesIRN), hold on
yline(min(dataDIP.lossesIRN))
figure, plot(dataDIP.lossesPEN)


%% draw ROIs
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(squeeze(dataDIP.T1Iter(20,:,:))));
axis image;
rois = get_rois;

%%

for i=1:size(dataDIP.T2Iter,1)
    for v=1:size(rois,2)
        [ix,iy] = find(rois{1,v}==1);
        t2iter{i,v}=reshape(squeeze(dataDIP.T2Iter(i,ix,iy)),[size(ix,1)*size(iy,1),1]);
        t2iter_m(i,v) = mean(t2iter{i,v});
        t2iter_std(i,v) = std(t2iter{i,v});
        
        t1iter{i,v}=reshape(squeeze(dataDIP.T1Iter(i,ix,iy)),[size(ix,1)*size(iy,1),1]);
        t1iter_m(i,v) = mean(t1iter{i,v});
        t1iter_std(i,v) = std(t1iter{i,v});


    end
end


%% import NIST data
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';
%%
iterstep=20;


fig=figure; 
fig.Position = [0 0 1000 1000];
t = tiledlayout(5, 4, "TileSpacing", "tight");
for v=1:size(rois,2)
    nexttile;
    errorbar([1:1:size(dataDIP.T2Iter,1)]*iterstep, t2iter_m(:,v), t2iter_std(:,v), 'o-'),
    hold on, yline(NISTplot(v,2),'r--', 'LineWidth',3)
    xlim([100,size(dataDIP.T2Iter,1)*iterstep])
end




fig=figure; 
fig.Position = [0 0 1000 1000];
t = tiledlayout(5, 4, "TileSpacing", "tight");
for v=1:size(rois,2)
    nexttile;
    errorbar([1:1:size(dataDIP.T1Iter,1)]*iterstep, t1iter_m(:,v), t1iter_std(:,v), 'o-'),
    hold on, yline(NISTplot(v,1),'r--', 'LineWidth',3)
    xlim([100,size(dataDIP.T1Iter,1)*iterstep])
end
