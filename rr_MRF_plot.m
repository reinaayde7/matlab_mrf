% rr MRF analysis
%considering one slice only!
clear all, close all, clc

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);
%%
MRFdir = ['E:\scanner_data\twix_data\241217\'];
MRFname = { 
    'MRF_fit_18-Dec-2024_meas_MID00026_BRAIN_s11_SVDdict_Nex1000'};

for i = 1: length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,i) = data{i}.t1map;
    t2map(:,:,i) = data{i}.t2map;
    m0map(:,:,i) = data{i}.m0map;
    if isfield(data{i}, 'b0map')
        b0map(:,:,i) = data{i}.b0map;
    else
        b0map(:,:,i) = zeros(size(data{i}.t1map));
    end

    % t1map(:,:,i) = data{i}.t1map3d(:,:,end);
    % t2map(:,:,i) = data{i}.t2map3d(:,:,end);
    % m0map(:,:,i) = data{i}.m0map3d(:,:,end);
    % b0map(:,:,i) = data{i}.b0map3d(:,:,end);
end

%% DIP results import

MRFDIPname{1} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\FISP\BRAIN11_241218_f300_ISO1mm_Nex1000_TE2_8_TR17_FAbody15\FISP_DIP_Drop10_Epochs200';
MRFDIPname{2} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\TRUEFISP\BRAIN11_241218_f300_ISO1mm_Nex1000_TE2_8_TR14_FAbody10\bSSFP_DIP_Drop5_Epochs200';
% MRFDIPname{2} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\FISP\v01_brain_s11_241018_t600\FISP_DIP_Drop5_Epochs200';
% MRFDIPname{3} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\TRUEFISP\v01_brain_s11_241018\bSSFP_DIP_Drop5_Epochs200';
% MRFDIPname{4} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\TRUEFISP\v01_brain_s11_241018_t600\bSSFP_DIP_Drop5_Epochs200';

for j = 1: length(MRFDIPname)
    dataDIP{j} = load([MRFDIPname{j} '.mat']);
    t1map(:,:,i+1) = squeeze(dataDIP{j}.T1Iter(end-2,:,:));
    t2map(:,:,i+1) = squeeze(dataDIP{j}.T2Iter(end-2,:,:));
    m0map(:,:,i+1) = squeeze(dataDIP{j}.M0Iter(end-2,:,:));
    if isfield(dataDIP{j}, 'B0Iter')
        b0map(:,:,i+1) = squeeze(dataDIP{j}.B0Iter(end,:,:));
    else
        b0map(:,:,i+1) = zeros(size(data{j}.t1map));
    end
    i=i+1;
end


%%

% full brain v01
s_vert = 60:240;
s_hor = 80:300;


% zoom in brain v01 GM v WM frontal
% s_vert = 210:270;
% s_hor = 270:330;

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 1500];
climst2 = [0, 300];
climstB0 = [-100 100];

for i = 1:length(MRFname)+length(MRFDIPname)
    f=figure('color','white');
    f.Position = [300 300 2400 420];
    ax1 = subplot(141); imagesc(rot90(t1map(s_vert,s_hor,i)), climst1), colormap(ax1, T1colormap), colorbar, axis off
    ax2 = subplot(142); imagesc(rot90(t2map(s_vert,s_hor,i)), climst2), colormap(ax2,T2colormap), colorbar, axis off
    mm = rot90(m0map(s_vert,s_hor,i));
    ax3 = subplot(143); imagesc(abs(mm)/max(abs(mm(:)))), colormap(ax3,gray), colorbar, axis off
    ax4 = subplot(144); imagesc(rot90(b0map(s_vert,s_hor,i)),climstB0), colormap(ax4,gray), colorbar, axis off
end

%%
