
% % phantom 
% folderdata = 'E:\scanner_data\twix_data\250203\';
% case88 = load([folderdata '3D_MRF_fit_04-Feb-2025_meas_MID00188_brain_FISP_FA15_nex600_image_combined_3Dfit_SVD_tilted.mat']);
% case58 = load([folderdata '3D_MRF_fit_04-Feb-2025_meas_MID00188_brain_FISP_FA15_nex600_image_combined_3Dfit_SVD_tilted_mockPF5_8.mat']);
% case68 = load([folderdata '3D_MRF_fit_04-Feb-2025_meas_MID00188_brain_FISP_FA15_nex600_image_combined_3Dfit_SVD_tilted_mockPF6_8.mat']);
% case78 = load([folderdata '3D_MRF_fit_05-Feb-2025_meas_MID00188_brain_FISP_FA15_nex600_image_combined_3Dfit_SVD_tilted_mockPF7_8.mat']);

% NIST
folderdata = 'E:\scanner_data\twix_data\250131\';
case88 = load([folderdata '3D_MRF_fit_03-Feb-2025_meas_MID00029_NIST_T2L_FISP_FA15_nex600_image_combined_3Dfit_SVD.mat']);
case58 = load([folderdata '3D_MRF_fit_03-Feb-2025_meas_MID00029_NIST_T2L_FISP_FA15_nex600_image_combined_3Dfit_SVD_PF58_mock.mat']);
case68 = load([folderdata '3D_MRF_fit_05-Feb-2025_meas_MID00029_NIST_T2L_FISP_FA15_nex600_image_combined_3Dfit_SVD_PF68_mock.mat']);
case78 = load([folderdata '3D_MRF_fit_05-Feb-2025_meas_MID00029_NIST_T2L_FISP_FA15_nex600_image_combined_3Dfit_SVD_PF78_mock.mat']);

%% T1
diff85 = squeeze(case88.t1map3d)-squeeze(case58.t1map3d);
diff86 = squeeze(case88.t1map3d)-squeeze(case68.t1map3d);
diff87 = squeeze(case88.t1map3d)-squeeze(case78.t1map3d);

line1 = cat(2,squeeze(case88.t1map3d),squeeze(case58.t1map3d), diff85);
line2 = cat(2,squeeze(case88.t1map3d),squeeze(case68.t1map3d), diff86);
line3 = cat(2,squeeze(case88.t1map3d),squeeze(case78.t1map3d), diff87);
group = cat(1,line1, line2, line3);

toSliceViewer = group;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%% T2
diff85 = squeeze(case88.t2map3d)-squeeze(case58.t2map3d);
diff86 = squeeze(case88.t2map3d)-squeeze(case68.t2map3d);
diff87 = squeeze(case88.t2map3d)-squeeze(case78.t2map3d);

line1 = cat(2,squeeze(case88.t2map3d),squeeze(case58.t2map3d), diff85);
line2 = cat(2,squeeze(case88.t2map3d),squeeze(case68.t2map3d), diff86);
line3 = cat(2,squeeze(case88.t2map3d),squeeze(case78.t2map3d), diff87);
group = cat(1,line1, line2, line3);
toSliceViewer = group;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

%% M0
case88.m0map = squeeze(abs(case88.m0map3d));
case58.m0map = squeeze(abs(case58.m0map3d));
case68.m0map = squeeze(abs(case68.m0map3d));
case78.m0map = squeeze(abs(case78.m0map3d));

diff85 = abs(squeeze(case88.m0map)-squeeze(case58.m0map));
diff86 = abs(squeeze(case88.m0map)-squeeze(case68.m0map));
diff87 = abs(squeeze(case88.m0map)-squeeze(case78.m0map));

line1 = cat(2,case88.m0map,case58.m0map, diff85);
line2 = cat(2,case88.m0map,case68.m0map, diff86);
line3 = cat(2,case88.m0map,case78.m0map, diff87);
group = cat(1,line1, line2, line3);
toSliceViewer = group;
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

%% import NIST data
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])

%% draw ROIs
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
sidx = 25;
figure;
imagesc(squeeze(case88.t1map3d(:,:,1,sidx)));
axis image;
rois = get_rois;

%% quantitative plots
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';

   
eval{1} = eval_vs_NIST(NISTplot, squeeze(case88.t1map3d(:,:,1,sidx)), squeeze(case88.t2map3d(:,:,1,sidx)), rois, squeeze(case88.m0map3d(:,:,1,sidx)));
eval{2} = eval_vs_NIST(NISTplot, squeeze(case58.t1map3d(:,:,1,sidx)), squeeze(case58.t2map3d(:,:,1,sidx)), rois, squeeze(case58.m0map3d(:,:,1,sidx)));
eval{2} = eval_vs_NIST(NISTplot, squeeze(case68.t1map3d(:,:,1,sidx)), squeeze(case68.t2map3d(:,:,1,sidx)), rois, squeeze(case68.m0map3d(:,:,1,sidx)));
eval{2} = eval_vs_NIST(NISTplot, squeeze(case78.t1map3d(:,:,1,sidx)), squeeze(case78.t2map3d(:,:,1,sidx)), rois, squeeze(case78.m0map3d(:,:,1,sidx)));

%%

difft1_1 = eval{1}.means_t1 - eval{2}.means_t1;
difft2_1 = eval{1}.means_t2 - eval{2}.means_t2;


figure,
subplot(221)
plot(eval{1}.means_t1, '-o'), hold on
plot(eval{2}.means_t1, '-o'), hold on
subplot(222)
plot(eval{1}.means_t2, '-o'), hold on
plot(eval{2}.means_t2, '-o'), hold on

subplot(223)
plot(difft1_1, '-o'), hold on
subplot(224)
plot(difft2_1, '-o'), hold on



%%
clearvars log
% evaluation excluding short T1s and T2s
th = 14;


% [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %
[mean_est_t1(:,1), mean_est_t2(:,1)] = REGRNIST(NISTplot, squeeze(case88.t1map3d(:,:,1,sidx)), squeeze(case88.t2map3d(:,:,1,sidx)), rois, [1,th]); %
[mean_est_t1(:,2), mean_est_t2(:,2)] = REGRNIST(NISTplot, squeeze(case58.t1map3d(:,:,1,sidx)), squeeze(case58.t2map3d(:,:,1,sidx)), rois, [1,th]); %
[mean_est_t1(:,3), mean_est_t2(:,3)] = REGRNIST(NISTplot, squeeze(case68.t1map3d(:,:,1,sidx)), squeeze(case68.t2map3d(:,:,1,sidx)), rois, [1,th]); %
[mean_est_t1(:,4), mean_est_t2(:,4)] = REGRNIST(NISTplot, squeeze(case78.t1map3d(:,:,1,sidx)), squeeze(case78.t2map3d(:,:,1,sidx)), rois, [1,th]); %

%regression
[coeffs_t1(1,:), R2_t1(1,:)]=regression_fit_values(log(NISTplot(1:th,1)), log(mean_est_t1(1:th,1)));
[coeffs_t2(2,:), R2_t2(2,:)]=regression_fit_values(log(NISTplot(1:th,2)), log(mean_est_t2(1:th,2)));

R2_t1 = R2_t1';
R2_t2 = R2_t2';

%%


err_tab_t1(:,1) = eval{1}.err_t1; err_tab_t1(:,2) = eval{2}.err_t1;
err_tab_t2(:,1) = eval{1}.err_t2; err_tab_t2(:,2) = eval{2}.err_t2;

std_tab_t1(:,1) = eval{1}.std_t1; std_tab_t1(:,2) = eval{2}.std_t1;
std_tab_t2(:,1) = eval{1}.std_t2; std_tab_t2(:,2) = eval{2}.std_t2;

cov_t1 = std_tab_t1./mean_est_t1;
cov_t2 = std_tab_t2./mean_est_t2;

rel_error_tab_t1 = err_tab_t1./repmat(NISTplot(:,1),[1,2])*100;
rel_error_tab_t2 = err_tab_t2./repmat(NISTplot(:,2),[1,2])*100;

mean_rel_error_t1_th = mean(rel_error_tab_t1(1:th,:)); %excluding some vials
mean_rel_error_t2_th = mean(rel_error_tab_t2(1:th,:));

clearvars err_tab_t2 err_tab_t1 std_tab_t1 std_tab_t2

%%
figure,
violinplot(cov_t2);

figure,
boxplot(cov_t2);

