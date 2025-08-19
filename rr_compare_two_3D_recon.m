addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500]; %1500/3000
climst2 = [0, 250]; %300/1200

folderdata = 'E:\scanner_data\twix_data\250228\';
case1 = load([folderdata '3D_MRF_fit_05-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD.mat']);
case2 = load([folderdata '3D_MRF_fit_05-Mar-2025_meas_MID00025_brain_HV4_FISP_FA15_nex600_w1500_3Dfit_SVD.mat']);

%% T1
toSliceViewer = cat(2,squeeze(case1.t1map3d),squeeze(case2.t1map3d),squeeze(case1.t1map3d)-squeeze(case2.t1map3d));
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%% T2
toSliceViewer = cat(2,squeeze(case1.t2map3d),squeeze(case2.t2map3d),squeeze(case1.t2map3d)-squeeze(case2.t2map3d));
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

%% M0
case1.m0map = squeeze(abs(case1.m0map3d));
case2.m0map = squeeze(abs(case2.m0map3d));

toSliceViewer = cat(2,squeeze(case1.m0map),squeeze(case2.m0map),abs(squeeze(case1.m0map)-squeeze(case2.m0map)));
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 0.2], 'ScaleFactors', [1,1,1]), colorbar

%% import NIST data
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])

%% draw ROIs
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
sidx = 23;
figure;
imagesc(squeeze(case1.t1map3d(:,:,1,sidx)));
axis image;
rois = get_rois;

%% quantitative plots
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';

   
eval{1} = eval_vs_NIST(NISTplot, squeeze(case1.t1map3d(:,:,1,sidx)), squeeze(case1.t2map3d(:,:,1,sidx)), rois, squeeze(case1.m0map3d(:,:,1,sidx)));
eval{2} = eval_vs_NIST(NISTplot, squeeze(case2.t1map3d(:,:,1,sidx)), squeeze(case2.t2map3d(:,:,1,sidx)), rois, squeeze(case2.m0map3d(:,:,1,sidx)));

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
[mean_est_t1(:,1), mean_est_t2(:,1)] = REGRNIST(NISTplot, squeeze(case1.t1map3d(:,:,1,sidx)), squeeze(case1.t2map3d(:,:,1,sidx)), rois, [1,th]); %
[mean_est_t1(:,2), mean_est_t2(:,2)] = REGRNIST(NISTplot, squeeze(case2.t1map3d(:,:,1,sidx)), squeeze(case2.t2map3d(:,:,1,sidx)), rois, [1,th]); %

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



%%
%%

for r =1:length(rois)
        [ix,iy] = find(rois{1,r}==1);
        fp=reshape(squeeze(case1.t1map3d(ix,iy,sidx)),[size(ix,1)*size(iy,1),1]);
        avgV(r,1) = mean(fp);
        stdV(r,1) = std(fp);

        fp=reshape(squeeze(case2.t1map3d(ix,iy,sidx)),[size(ix,1)*size(iy,1),1]);
        avgV(r,2) = mean(fp);
        stdV(r,2) = std(fp);

        fp=reshape(squeeze(case1.t2map3d(ix,iy,sidx)),[size(ix,1)*size(iy,1),1]);
        avgV(r,3) = mean(fp);
        stdV(r,3) = std(fp);

        fp=reshape(squeeze(case2.t2map3d(ix,iy,sidx)),[size(ix,1)*size(iy,1),1]);
        avgV(r,4) = mean(fp);
        stdV(r,4) = std(fp);
end

%%
% to do: add manual skull stripping.
%masking
% T2
for slice = 1:48
    [idx_x, idx_y] = ind2sub(size(case1.t1map3d(:,:,1,slice)), find(case1.t1map3d(:,:,1,slice) >= 600 & case1.t1map3d(:,:,1,slice) <= 1300)); %GM
    [idx_x, idx_y] = ind2sub(size(case1.t1map3d(:,:,1,slice)), find(case1.t1map3d(:,:,1,slice) >= 0 & case1.t1map3d(:,:,1,slice) <= 600)); %WM
    [idx_x, idx_y] = ind2sub(size(case1.t1map3d(:,:,1,slice)), find(case1.t1map3d(:,:,1,slice) >= 1300)); %CSF
    m = zeros(size(case1.t1map3d,1), size(case1.t1map3d,2)); 

    % Assign 1 to the corresponding (y, x) positions
    m(sub2ind(size(mask), idx_x, idx_y)) = 1;
    mask(:,:,slice) = m;
end

%%
toSliceViewer = cat(2,squeeze(case1.t1map3d).*mask,squeeze(case2.t1map3d).*mask,(squeeze(case1.t1map3d)-squeeze(case2.t1map3d)).*mask);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar