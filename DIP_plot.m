% rr MRF analysis
%considering one slice only!
clear all, close all, clc

%%
MRFdir = ['E:\scanner_data\twix_data\241011\'];
MRFname = { 'Rudy_NIST_2024_1024_DIP_Drop5_Epochs300';
    };

% MRFname = { 'Rudy_VIVO_2024_1024_DIP_Drop5_Epochs300';
%     };
%%
data = load([MRFdir MRFname{1} '.mat']);
%%

t1map = squeeze(data.T1Iter(end,:,:));
t2map = squeeze(data.T2Iter(end,:,:));
m0map = squeeze(data.M0Iter(end,:,:));
b0map = squeeze(data.DFIter(end,:,:));

%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 3000];
climst2 = [0, 300]; %uplimit:1200 for NIST T2 at 0.55T 
climsb0 = [-40,40];

%%
figure,
ax1 = subplot(131); imagesc(t1map, climst1), colormap(ax1, T1colormap), colorbar
ax2 = subplot(132); imagesc(t2map, climst2), colormap(ax2,T2colormap), colorbar
ax3 = subplot(133); imagesc(abs(m0map)/max(abs(m0map(:)))), colormap(ax3,gray), colorbar

figure, imagesc(b0map,climsb0), colormap(gray), colorbar

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
for i = 1: length(MRFname)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois);
    eval{i} = eval_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois, b0map(:,:,i));
end
%%
for i = 1: length(MRFname)
    [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %
end

%%
err_tab_t1 = zeros(14,size(data,2));
err_tab_t2 = zeros(14,size(data,2));
std_tab_t1 = zeros(14,size(data,2));
std_tab_t2 = zeros(14,size(data,2));
for i=1:size(data,2)
    err_tab_t1(:,i) = eval{i}.err_t1;
    err_tab_t2(:,i) = eval{i}.err_t2;

    std_tab_t1(:,i) = eval{i}.std_t1;
    std_tab_t2(:,i) = eval{i}.std_t2;
end

%%
figure, h = heatmap(abs(err_tab_t2), 'Colormap', turbo);
h.YData = NISTplot(:,2);
h.Title = 'T2 error - heatmap';

figure, h = heatmap(abs(err_tab_t1), 'Colormap', turbo);
h.YData = NISTplot(:,1);
h.Title = 'T1 error - heatmap';

figure, h = heatmap(abs(std_tab_t2), 'Colormap', turbo);
h.YData = NISTplot(:,2);
h.Title = 'T2 std - heatmap';

figure, h = heatmap(abs(std_tab_t1), 'Colormap', turbo);
h.YData = NISTplot(:,1);
h.Title = 'T1 std - heatmap';
