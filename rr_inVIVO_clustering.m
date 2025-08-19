
clear all, close all, clc

%HV 1
MRFdir = ['E:\scanner_data\twix_data\250203\'];
MRFname = { '3D_MRF_fit_07-Apr-2025_meas_MID00188_brain_HV1_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
    '3D_MRF_DIP_12-Feb-2025_MID00188_DIP_3Dvolume' %DIP bSSFP 1000
    };

% %HV 2
% MRFdir = ['E:\scanner_data\twix_data\250221\'];
% MRFname = { '3D_MRF_fit_31-Mar-2025_meas_MID00547_brain_HV2_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
%     '3D_MRF_DIP_03-Mar-2025_meas_MID00547_brain_HV2_FISP_FA15_nex600_DIP_3Dvolume_HV2' %DIP bSSFP 1000
%     };

%HV 3
% MRFdir = ['E:\scanner_data\twix_data\250221\'];
% MRFname = { '3D_MRF_fit_31-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
%     '3D_MRF_DIP_03-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_DIP_3Dvolume_HV3' %DIP bSSFP 1000
%     };

% %HV 4
% MRFdir = ['E:\scanner_data\twix_data\250228\'];
% MRFname = { '3D_MRF_fit_07-Apr-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
%     '3D_MRF_DIP_07-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD'
%     };

%HV 5
% MRFdir = ['E:\scanner_data\twix_data\250305\'];
% MRFname = { '3D_MRF_fit_07-Apr-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
%     '3D_MRF_DIP_07-Mar-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD'
%     };
% 
%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 250]; %uplimit:1200 for NIST T2 at 0.55T 

%LR
data{1} = load([MRFdir MRFname{1} '.mat']);
t1map(:,:,:,1) = rot90(squeeze(data{1}.t1map3d));
t2map(:,:,:,1) = rot90(squeeze(data{1}.t2map3d));
m0map(:,:,:,1) = rot90(squeeze(data{1}.m0map3d));
b0map(:,:,:,1) = rot90(squeeze(data{1}.b0map3d));

%SLLR
t1map(:,:,:,2) = rot90(squeeze(data{1}.t1map3dLLR));
t2map(:,:,:,2) = rot90(squeeze(data{1}.t2map3dLLR));
m0map(:,:,:,2) = rot90(squeeze(data{1}.m0map3dLLR));
b0map(:,:,:,2) = rot90(squeeze(data{1}.b0map3dLLR));

%DIP
data{2} = load([MRFdir MRFname{2} '.mat']);
t1map(:,:,:,3) = squeeze(data{2}.t1map3d);
t2map(:,:,:,3) = squeeze(data{2}.t2map3d);
m0map(:,:,:,3) = squeeze(data{2}.m0map3d);
b0map(:,:,:,3) = rot90(squeeze(data{2}.b0map3d));

toSliceViewer = squeeze(t1map(:,:,:,1));


figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar


%% automatic selection of WM/GM and CSF volumes
% select slices superior to ventricles
% manual skull stripping
% apply 3D k-means clustering with 5 clusters
% identify clusters for tissue prop by visual inspection
% morphological erosion applied to WM to eliminate partial voluming
% evaluate T1 and T2 in masked regions

olds = load('E:\POSTDOC_UoM\10_manuscripts\DIP_brainMRF_055T\data_HV_cluster\HV1_ROI_stats_clusterLR.mat');
%%
slices_for_analysis = olds.slices_for_analysis; 
skull_mask = olds.skull_mask;
brain_t1maps = t1map(:,:,slices_for_analysis,:).*skull_mask;
toSliceViewer = squeeze(brain_t1maps(:,:,:,2));

figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%% freehand masking volunteer specific - skull stripping
% select 10 slices superior to the ventricles
slices_for_analysis = [25:35]; 

codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
for n = 1:numel(slices_for_analysis)
    slice = slices_for_analysis(n);
    figure;
    imagesc(abs(squeeze(t1map(:,:,slice,2))), climst1);
    colormap(T1colormap)
    axis image;
    rois{n} = get_rois;
end

for n = 1:numel(slices_for_analysis)
    skull_mask(:,:,n) = rois{n}{1}; 
end

brain_t1maps = t1map(:,:,slices_for_analysis,:).*skull_mask;

toSliceViewer = squeeze(brain_t1maps(:,:,:,1));

figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%% kmeans
T1_matrix = brain_t1maps(:,:,:,2); %cluster on the LR model!
se = strel('sphere', 1); % 3D spherical structuring element with radius 2 voxels

% T1_matrix = brain_t1maps(:,:,:,2); %cluster on the DIP model!
% se = strel('sphere', 2); % 3D spherical structuring element with radius 2 voxels


% Assume T1_matrix is your 3D matrix of T1 values
k = 5; % Number of clusters

% Get the size of the original matrix
[xDim, yDim, zDim] = size(T1_matrix);

% Reshape 3D matrix into a 2D matrix where each row is a voxel
data = T1_matrix(:); % Convert to a column vector

% Run k-means clustering
[idx, C] = kmeans(data, k, 'MaxIter', 1000, 'Replicates', 5);

% Reshape idx back to 3D to match the original matrix
cluster_labels = reshape(idx, [xDim, yDim, zDim]);

% Create 5 binary masks (one for each cluster)
cluster_masks = cell(k,1); % Store each cluster as a separate 3D matrix
for i = 1:k
    cluster_masks{i} = (cluster_labels == i); % Binary mask for cluster i
end

toSliceViewer = cat(2, cluster_masks{1},cluster_masks{2},cluster_masks{3},cluster_masks{4},cluster_masks{5});
figure,
sliceViewer(toSliceViewer), colorbar

%%
% Select the cluster mask you want to erode (e.g., Cluster 1) -> WM cluster
cluster_to_erode = cluster_masks{5}; 

% Define a 3D structuring element (a spherical or cubic kernel)
%done above

% Apply erosion
eroded_cluster = imerode(cluster_to_erode, se);

% Visualize the middle slice before and after erosion
zSlice = round(size(cluster_to_erode, 3) / 2);
figure;
subplot(1,2,1);
imagesc(cluster_to_erode(:,:,zSlice)); 
colormap gray; axis equal; title('Original Cluster Mask');

subplot(1,2,2);
imagesc(eroded_cluster(:,:,zSlice)); 
colormap gray; axis equal; title('Eroded Cluster Mask');


%% sanity check
% wm_cluster = eroded_cluster;
% gm_cluster = cluster_masks{1};% +cluster_masks{3};
% csf_cluster = cluster_masks{3};

wm_cluster = olds.wm_cluster;
gm_cluster = olds.gm_cluster;
csf_cluster = olds.csf_cluster;
cluster_masks = olds.cluster_masks;

toSliceViewer = cat(2, wm_cluster,gm_cluster,csf_cluster);
figure,
sliceViewer(toSliceViewer), colorbar

%%
T1_WM = brain_t1maps(:,:,:,1).*wm_cluster;
T1_GM = brain_t1maps(:,:,:,1).*gm_cluster;
T1_CSF = brain_t1maps(:,:,:,1).*csf_cluster;

T1_WM_llr = brain_t1maps(:,:,:,2).*wm_cluster;
T1_GM_llr = brain_t1maps(:,:,:,2).*gm_cluster;
T1_CSF_llr = brain_t1maps(:,:,:,2).*csf_cluster;

T1_WM_dip = brain_t1maps(:,:,:,3).*wm_cluster;
T1_GM_dip = brain_t1maps(:,:,:,3).*gm_cluster;
T1_CSF_dip = brain_t1maps(:,:,:,3).*csf_cluster;

brain_t2maps = t2map(:,:,slices_for_analysis,:).*skull_mask;

T2_WM = brain_t2maps(:,:,:,1).*wm_cluster;
T2_GM = brain_t2maps(:,:,:,1).*gm_cluster;
T2_CSF = brain_t2maps(:,:,:,1).*csf_cluster;

T2_WM_llr = brain_t2maps(:,:,:,2).*wm_cluster;
T2_GM_llr = brain_t2maps(:,:,:,2).*gm_cluster;
T2_CSF_llr = brain_t2maps(:,:,:,2).*csf_cluster;

T2_WM_dip = brain_t2maps(:,:,:,3).*wm_cluster;
T2_GM_dip = brain_t2maps(:,:,:,3).*gm_cluster;
T2_CSF_dip = brain_t2maps(:,:,:,3).*csf_cluster;


%%
toSliceViewer = cat(2, T1_WM,T1_WM_llr,T1_WM_dip,T1_GM,T1_GM_llr,T1_GM_dip);
figure,
sliceViewer(toSliceViewer), colorbar
%%
T1v_WM = T1_WM(:); T1v_WM_llr = T1_WM_llr(:); T1v_WM_dip = T1_WM_dip(:);
T1v_WM = T1v_WM(T1v_WM>20); T1v_WM_llr = T1v_WM_llr(T1v_WM_llr>20);  T1v_WM_dip = T1v_WM_dip(T1v_WM_dip>20);

T2v_WM = T2_WM(:);T2v_WM_llr = T2_WM_llr(:); T2v_WM_dip = T2_WM_dip(:);
T2v_WM = T2v_WM(T2v_WM>5); T2v_WM_llr = T2v_WM_llr(T2v_WM_llr>5); T2v_WM_dip = T2v_WM_dip(T2v_WM_dip>5);

figure, 
subplot(121), histogram(T1v_WM,100), hold on, histogram(T1v_WM_llr,100), hold on, histogram(T1v_WM_dip,100)
subplot(122), histogram(T2v_WM,100), hold on, histogram(T2v_WM_llr,100), hold on, histogram(T2v_WM_dip,100)

%%
T1v_GM = T1_GM(:); T1v_GM_llr = T1_GM_llr(:); T1v_GM_dip = T1_GM_dip(:);
T1v_GM = T1v_GM(T1v_GM>100); T1v_GM_llr = T1v_GM_llr(T1v_GM_llr>100); T1v_GM_dip = T1v_GM_dip(T1v_GM_dip>100);

T2v_GM = T2_GM(:); T2v_GM_llr = T2_GM_llr(:); T2v_GM_dip = T2_GM_dip(:);
T2v_GM = T2v_GM(T2v_GM>20); T2v_GM_llr = T2v_GM_llr(T2v_GM_llr>20); T2v_GM_dip = T2v_GM_dip(T2v_GM_dip>20);

figure, 
subplot(121), histogram(T1v_GM,100),hold on, histogram(T1v_GM_llr,100),hold on, histogram(T1v_GM_dip,100)
subplot(122), histogram(T2v_GM,100),hold on, histogram(T2v_GM_llr,100),hold on, histogram(T2v_GM_dip,100)

%%
T1v_CSF = T1_CSF(:); T1v_CSF_llr = T1_CSF_llr(:); T1v_CSF_dip = T1_CSF_dip(:);
T1v_CSF = T1v_CSF(T1v_CSF>100); T1v_CSF_llr = T1v_CSF_llr(T1v_CSF_llr>100); T1v_CSF_dip = T1v_CSF_dip(T1v_CSF_dip>100);

T2v_CSF = T2_CSF(:); T2v_CSF_llr = T2_CSF_llr(:); T2v_CSF_dip = T2_CSF_dip(:);
T2v_CSF = T2v_CSF(T2v_CSF>20); T2v_CSF_llr = T2v_CSF_llr(T2v_CSF_llr>20); T2v_CSF_dip = T2v_CSF_dip(T2v_CSF_dip>20);

figure, 
subplot(121), histogram(T1v_CSF,100), hold on, histogram(T1v_CSF_llr,100), hold on, histogram(T1v_CSF_dip,100)
subplot(122), histogram(T2v_CSF,100), hold on, histogram(T2v_CSF_llr,100), hold on, histogram(T2v_CSF_dip,100)


%%
% Compute mean and standard deviation
mean_T1 = median(T1v_WM);mean_T1llr= median(T1v_WM_llr);mean_T1dip = median(T1v_WM_dip);
std_T1 = std(T1v_WM);std_T1llr = std(T1v_WM_llr);std_T1dip = std(T1v_WM_dip);
% Display with 2 decimal precision
fprintf('LR - T1_WM: %.2f pm %.2f \n', mean_T1, std_T1);
fprintf('SLLR - T1_WM: %.2f pm %.2f \n', mean_T1llr, std_T1llr);
fprintf('DIP - T1_WM: %.2f pm %.2f \n', mean_T1dip, std_T1dip);


mean_T1 = median(T1v_GM); mean_T1llr = median(T1v_GM_llr); mean_T1dip = median(T1v_GM_dip);
std_T1 = std(T1v_GM);std_T1llr = std(T1v_GM_llr);std_T1dip = std(T1v_GM_dip);
fprintf('LR - T1_GM: %.2f pm %.2f \n', mean_T1, std_T1);
fprintf('SLLR - T1_GM: %.2f pm %.2f \n', mean_T1llr, std_T1llr);
fprintf('DIP - T1_GM: %.2f pm %.2f \n', mean_T1dip, std_T1dip);

mean_T1 = median(T1v_CSF);mean_T1llr = median(T1v_CSF_llr);mean_T1dip = median(T1v_CSF_dip);
std_T1 = std(T1v_CSF);std_T1llr = std(T1v_CSF_llr);std_T1dip = std(T1v_CSF_dip);
fprintf('LR - T1_CSF: %.2f pm %.2f \n', mean_T1, std_T1);
fprintf('SLLR - T1_CSF: %.2f pm %.2f \n', mean_T1llr, std_T1llr);
fprintf('DIP - T1_CSF: %.2f pm %.2f \n', mean_T1dip, std_T1dip);

fprintf('----------------------------------------- \n');

% Compute mean and standard deviation
mean_T2 = median(T2v_WM);mean_T2llr = median(T2v_WM_llr);mean_T2dip = median(T2v_WM_dip);
std_T2 = std(T2v_WM);std_T2llr = std(T2v_WM_llr);std_T2dip = std(T2v_WM_dip);
% Display with 2 decimal precision
fprintf('LR - T2_WM: %.2f pm %.2f \n', mean_T2, std_T2);
fprintf('SLLR - T2_WM: %.2f pm %.2f \n', mean_T2llr, std_T2llr);
fprintf('DIP - T2_WM: %.2f pm %.2f \n', mean_T2dip, std_T2dip);

mean_T2 = median(T2v_GM); mean_T2llr = median(T2v_GM_llr); mean_T2dip = median(T2v_GM_dip);
std_T2 = std(T2v_GM);std_T2llr = std(T2v_GM_llr);std_T2dip = std(T2v_GM_dip);
fprintf('LR - T2_GM: %.2f pm %.2f \n', mean_T2, std_T2);
fprintf('SLLR - T2_GM: %.2f pm %.2f \n', mean_T2llr, std_T2llr);
fprintf('DIP - T2_GM: %.2f pm %.2f \n', mean_T2dip, std_T2dip);

mean_T2 = median(T2v_CSF);mean_T2llr = median(T2v_CSF_llr);mean_T2dip = median(T2v_CSF_dip);
std_T2 = std(T2v_CSF);std_T2llr = std(T2v_CSF_llr);std_T2dip = std(T2v_CSF_dip);
fprintf('LR - T2_CSF: %.2f pm %.2f \n', mean_T2, std_T2);
fprintf('SLLR - T2_CSF: %.2f pm %.2f \n', mean_T2llr, std_T2llr);
fprintf('DIP - T2_CSF: %.2f pm %.2f \n', mean_T2dip, std_T2dip);

%% savings
% saveName = 'HV5_ROI_stats_clusterLR.mat';
saveName = 'HV1_ROI_stats_clusterLR.mat';
save([MRFdir saveName], 'T1_WM', 'T1_WM_llr', 'T1_WM_dip', 'T1_GM', 'T1_GM_llr', 'T1_GM_dip', 'T1_CSF', 'T1_CSF_llr', 'T1_CSF_dip', ...
    'T2_WM', 'T2_WM_llr', 'T2_WM_dip', 'T2_GM', 'T2_GM_llr', 'T2_GM_dip', 'T2_CSF', 'T2_CSF_llr', 'T2_CSF_dip', ...
    'cluster_masks', 'wm_cluster', 'gm_cluster' ,'csf_cluster', 'brain_t1maps', 'brain_t2maps','slices_for_analysis','skull_mask');
