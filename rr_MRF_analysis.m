% rr MRF analysis
%considering one slice only!
clear all, close all, clc

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);
%%
MRFdir = ['E:\scanner_data\twix_data\250328\'];
MRFname = { 
    % 'MRF_fit_28-Mar-2025_meas_MID00200_NIST_s20_SVDdict_Nex600_FA15_FISP_w1500_LR_LLR';
    % 'MRF_fit_28-Mar-2025_meas_MID00202_NIST_s20_SVDdict_Nex600_FAopt_FISP_w1500_LR_LLR';
    'MRF_fit_31-Mar-2025_meas_MID00201_NIST_s20_SVDdict_Nex600_FA15_pSSFP_w1500_LR_LLR';
    'MRF_fit_01-Apr-2025_meas_MID00201_NIST_s20_SVDdict_Nex600_FAopt_pSSFP_w1500_b0_0_LR_LLR';
    };

k=1;
for i = 1:length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,k) = data{i}.t1map;
    t2map(:,:,k) = data{i}.t2map;
    m0map(:,:,k) = data{i}.m0map;
    if isfield(data{i}, 'b0map')
        b0map(:,:,k) = data{i}.b0map;
    else
        b0map(:,:,k) = zeros(size(data{i}.t1map));
    end
    if isfield(data{i}, 't1mapLLR')
        k=k+1;
        t1map(:,:,k) = data{i}.t1mapLLR;
        t2map(:,:,k) = data{i}.t2mapLLR;
        m0map(:,:,k) = data{i}.m0mapLLR;
        if isfield(data{i}, 'b0map')
            b0map(:,:,k) = data{i}.b0mapLLR;
        else
            b0map(:,:,k) = zeros(size(data{i}.t1mapLLR));
        end
    end
    k=k+1;

end

%% DIP results import

% dictionary gen on the v001 - NIST: FISP 1000v600
% MRFDIPname{1} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\FISP\NIST_T2L_250328_s20_f300_Nex600_FAbody15_wait1500\FISP_DIP_Drop5_Epochs400';
% MRFDIPname{2} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\FISP\NIST_T2L_250328_s20_f300_Nex600_FAopt_wait1500\FISP_DIP_Drop5_Epochs400';
MRFDIPname{1} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\TRUEFISP\NIST_T2L_250328_s20_f300_Nex600_FA15_wait1500\bSSFP_DIP_Drop5_Epochs400';
MRFDIPname{2} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\TRUEFISP\NIST_T2L_250328_s20_f300_Nex600_FAopt_wait1500\bSSFP_DIP_Drop5_Epochs400';
%MRFDIPname = {};
%%
for j = 1: length(MRFDIPname)
    dataDIP{j} = load([MRFDIPname{j} '.mat']);
    t1map(:,:,k) = squeeze(dataDIP{j}.T1Iter(10,:,:));
    t2map(:,:,k) = squeeze(dataDIP{j}.T2Iter(10,:,:));
    m0map(:,:,k) = squeeze(dataDIP{j}.M0Iter(10,:,:));
    if isfield(dataDIP{j}, 'B0Iter')
        b0map(:,:,k) = squeeze(dataDIP{j}.B0Iter(10,:,:));
    else
        b0map(:,:,k) = zeros(size(t1map,1), size(t1map,2));
    end
    k=k+1;
end

%% import NIST data
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])

%% draw ROIs
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(squeeze(t1map(:,:,1))));
axis image;
rois = get_rois;
%% quantitative plots
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';
for i =1:size(t1map,3)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois);
    eval{i} = eval_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois, b0map(:,:,i));
end
%%
% evaluation excluding short T1s and T2s
th_t1 = [1,14];
th_t2 = [1,14];
for i = 1:size(t1map,3)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %
    % [mean_est_t1(:,i), mean_est_t2(:,i)] = REGRNIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois, [1,th]); %
    [mean_est_t1(:,i), mean_est_t2(:,i)] = REGRNIST(NISTplot, eval{i}, th_t1,th_t2); %
    
end
%%
labels{1}='LR-FA15';
labels{2}='DIP-FA15';
[eo1, eo3] = REGRNIST_1v1(NISTplot, eval{3}, eval{4}, [4,11], labels); %
% 

%%
labels{1}='LR';
labels{2}='SLLR';
labels{3}='DIP';
[eo1, eo2, eo5] = REGRNIST_1v1v1(NISTplot, eval{3}, eval{4}, eval{6}, [1,14], labels); %


%%
FISP_LR = eo1;
FISP_SLLR = eo2;
FISP_DIP = eo5;
% pSSFP_LR = eo2;
% pSSFP_DIP = eo4;

%%
savedir = 'E:\POSTDOC_UoM\10_manuscripts\DIP_brainMRF_055T\';
save([savedir 'eval_FISP_LR_SLLR_DIP_FA15.mat'], 'FISP_LR','FISP_DIP','FISP_SLLR')

%%
s_vert = 50:250;
s_hor = 80:270;

% t1mapN=t1map;
% t1mapN(:,:,2) = t1map(:,:,3);
% t1mapN(:,:,3) = t1map(:,:,2);
% t2mapN=t2map;
% t2mapN(:,:,2) = t2map(:,:,3);
% t2mapN(:,:,3) = t2map(:,:,2);
% m0mapN=m0map;
% m0mapN(:,:,2) = m0map(:,:,3);
% m0mapN(:,:,3) = m0map(:,:,2);
% b0mapN=b0map;
% b0mapN(:,:,2) = b0map(:,:,3);
% b0mapN(:,:,3) = b0map(:,:,2);
% 
% t1map = t1mapN; t2map = t2mapN; m0map=m0mapN; b0map=b0mapN;
%%
for i=1:size(m0map,3)
    cc = abs(squeeze(m0map(:,:,i)));
    m0(:,:,i) = cc/max(cc(:));
end
%% mask calculation
figure, histogram(m0(:,:,5))
th = 0.03;

M0s = m0(:,:,5);
cluster_mask = M0s>th;


figure,
imagesc(cluster_mask), colorbar


%%

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 3000];
climst2 = [0, 1500];
climstB0 = [-100 100];
fig=figure; 
fig.Position = [0 0 1600 800];
t = tiledlayout(3, 6, "TileSpacing", "none");

for i = 1:size(t1map,3)
    ax = nexttile;
    plot = t1map(:,:,i).*cluster_mask;
    imagesc(plot(s_vert,s_hor), climst1), colormap(ax, T1colormap), axis off 
end
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;
    

for i = 1:size(t1map,3)
    ax = nexttile;
     plot = t2map(:,:,i).*cluster_mask;
    imagesc(plot(s_vert,s_hor), climst2), colormap(ax, T2colormap), axis off 
end
cb1 = colorbar(ax);
cb1.Position(4) = cb1.Position(4)*0.9;
    
for i = 1:size(t1map,3)
    ax = nexttile;
    plot = m0(:,:,i).*cluster_mask;
    imagesc(plot(s_vert,s_hor), [0 1]), colormap(ax, gray), axis off 
end

cb2 = colorbar(ax);
cb2.Position(4) = cb2.Position(4)*0.9;

% for i = 1:length(MRFname)+length(MRFDIPname)
%     ax = nexttile;
%     imagesc(b0map(s_vert,s_hor,i),climstB0), colormap(ax, gray), axis off 
% end
% 
% % nexttile; axis off; nexttile; axis off;
% % ax = nexttile;
% % imagesc(b0map(s_vert,s_hor,3), climstB0), colormap(ax, gray), axis off 
% % ax = nexttile;
% % imagesc(b0map(s_vert,s_hor,4), climstB0), colormap(ax, gray), axis off
% cb3 = colorbar(ax);
% cb3.Position(4) = cb3.Position(4)*0.9;

%%
climstB0 = [-100 100];
fig=figure; 
fig.Position = [0 0 1600 300];
t = tiledlayout(1, 6, "TileSpacing", "none");

for i = 1:size(t1map,3)
    ax = nexttile;
    plot = t1map(:,:,i).*cluster_mask;
    imagesc(b0map(s_vert,s_hor,i),climstB0), colormap(ax, gray), axis off 
end
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;