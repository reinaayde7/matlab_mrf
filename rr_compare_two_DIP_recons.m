% rr MRF analysis
%considering one slice only!
clear all, close all, clc

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);

%% DIP results import

% dictionary gen on the v001 - NIST: FISP 1000v600
DIPfolder = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\FISP\';
MRFDIPname{1} = 'NIST_T2L_250225_s23_f300_Nex600_w3000_FAbody10\FISP_DIP_Drop5_Epochs500';
MRFDIPname{2} = 'NIST_T2L_250225_s23_f300_Nex600_w3000_FAbody15\FISP_DIP_Drop5_Epochs500';

% MRFDIPname = {};

for j = 1: length(MRFDIPname)
    dataDIP{j} = load([DIPfolder MRFDIPname{j} '.mat']);
    t1map(:,:,j) = squeeze(dataDIP{j}.T1Iter(10,:,:));
    t2map(:,:,j) = squeeze(dataDIP{j}.T2Iter(10,:,:));
    m0map(:,:,j) = squeeze(dataDIP{j}.M0Iter(10,:,:));
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
for i =1:length(MRFDIPname)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois);
    eval{i} = eval_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois);
end
%%
% evaluation excluding short T1s and T2s
th = 11;
for i = 1:length(MRFDIPname)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %
    [mean_est_t1(:,i), mean_est_t2(:,i)] = REGRNIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois, [1,th]); %

    
    %regression
    [coeffs_t1(i,:), R2_t1(i)]=regression_fit_values(log(NISTplot(1:th, 1)), log(mean_est_t1(1:th,i)));
    [coeffs_t2(i,:), R2_t2(i)]=regression_fit_values(log(NISTplot(1:th, 2)), log(mean_est_t2(1:th,i)));
end

R2_t1 = R2_t1';
R2_t2 = R2_t2';
%%
% err_tab_t1 = zeros(14,size(data,2));
% err_tab_t2 = zeros(14,size(data,2));
% std_tab_t1 = zeros(14,size(data,2));
% std_tab_t2 = zeros(14,size(data,2));
for i=1:length(MRFDIPname)
    err_tab_t1(:,i) = eval{i}.err_t1;
    err_tab_t2(:,i) = eval{i}.err_t2;

    std_tab_t1(:,i) = eval{i}.std_t1;
    std_tab_t2(:,i) = eval{i}.std_t2;
end

cov_t1 = std_tab_t1./mean_est_t1;
cov_t2 = std_tab_t2./mean_est_t2;

rel_error_tab_t1 = err_tab_t1./repmat(NISTplot(:,1),[1,length(MRFDIPname)])*100;
rel_error_tab_t2 = err_tab_t2./repmat(NISTplot(:,2),[1,length(MRFDIPname)])*100;

mean_rel_error_t1_th = mean(rel_error_tab_t1(1:th,:)); %excluding some vials
mean_rel_error_t2_th = mean(rel_error_tab_t2(1:th,:));

clearvars err_tab_t2 err_tab_t1 std_tab_t1 std_tab_t2
%% plots over iter
for meas = 1: length(MRFDIPname)
    for i=1:size(dataDIP{1}.T2Iter,1)
        for v=1:size(rois,2)
            [ix,iy] = find(rois{1,v}==1);
            t2iter{i,v,meas}=reshape(squeeze(dataDIP{meas}.T2Iter(i,ix,iy)),[size(ix,1)*size(iy,1),1]);
            t2iter_m(i,v,meas) = mean(t2iter{i,v,meas});
            t2iter_std(i,v,meas) = std(t2iter{i,v,meas});
            
            t1iter{i,v,meas}=reshape(squeeze(dataDIP{meas}.T1Iter(i,ix,iy)),[size(ix,1)*size(iy,1),1]);
            t1iter_m(i,v,meas) = mean(t1iter{i,v,meas});
            t1iter_std(i,v,meas) = std(t1iter{i,v,meas});
    
    
        end
    end
end

%%
iterstep=20;

for meas = 1: length(MRFDIPname)
    fig=figure; 
    fig.Position = [0 0 1000 1000];
    t = tiledlayout(5, 4, "TileSpacing", "tight");
    for v=1:size(rois,2)
        nexttile;
        errorbar([1:1:size(dataDIP{meas}.T2Iter,1)]*iterstep, t2iter_m(:,v,meas), t2iter_std(:,v,meas), 'o-'),
        hold on, yline(NISTplot(v,2),'r--', 'LineWidth',3)
        xlim([100,size(dataDIP{meas}.T2Iter,1)*iterstep])
    end
end


for meas = 1: length(MRFDIPname)
    fig=figure; 
    fig.Position = [0 0 1000 1000];
    t = tiledlayout(5, 4, "TileSpacing", "tight");
    for v=1:size(rois,2)
        nexttile;
        errorbar([1:1:size(dataDIP{meas}.T1Iter,1)]*iterstep, t1iter_m(:,v,meas), t1iter_std(:,v,meas), 'o-'),
        hold on, yline(NISTplot(v,1),'r--', 'LineWidth',3)
        xlim([100,size(dataDIP{meas}.T1Iter,1)*iterstep])
    end
end
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

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 3000];
climst2 = [0, 1500];
climstB0 = [-40 40];
fig=figure; 
fig.Position = [0 0 1600 1000];
t = tiledlayout(4, 7, "TileSpacing", "none");

for i = 1:length(MRFDIPname)
    ax = nexttile;
    imagesc(t1map(s_vert,s_hor,i), climst1), colormap(ax, T1colormap), axis off 
end
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;
    

for i = 1:length(MRFDIPname)
    ax = nexttile;
    imagesc(t2map(s_vert,s_hor,i), climst2), colormap(ax, T2colormap), axis off 
end
cb1 = colorbar(ax);
cb1.Position(4) = cb1.Position(4)*0.9;
    
for i = 1:length(MRFDIPname)
    ax = nexttile;
    mm = m0map(s_vert,s_hor,i);
    imagesc(abs(mm)/max(abs(mm(:))), [0 1]), colormap(ax, gray), axis off 
end

cb2 = colorbar(ax);
cb2.Position(4) = cb2.Position(4)*0.9;

%%
% ax2 = subplot(142); imagesc(t2map(s_vert,s_hor,i), climst2), colormap(ax2,T2colormap), colorbar, axis off
    % mm = m0map(s_vert,s_hor,i);
    % ax3 = subplot(143); imagesc(abs(mm)/max(abs(mm(:)))), colormap(ax3,gray), colorbar, axis off
    % ax4 = subplot(144); imagesc(b0map(s_vert,s_hor,i),climstB0), colormap(ax4,gray), colorbar, axis off