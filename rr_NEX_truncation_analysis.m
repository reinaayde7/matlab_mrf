% rr MRF analysis
%considering one slice only!
clear all, close all, clc

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);

%% DIP results import

%FISP
MRFDIPname{1} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\FISP\NIST_241011\FISP_DIP_Drop5_Epochs200'; %Nex 1000
MRFDIPname{2} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\FISP\NIST_241011_t600\FISP_DIP_Drop5_Epochs200'; %Nex600
MRFDIPname{3} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\FISP\NIST_241011_t300\FISP_DIP_Drop5_Epochs200'; %Nex300

%trueFISP
MRFDIPname{4} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\TRUEFISP\NIST_241011\bSSFP_DIP_Drop5_Epochs200'; %Nex 1000
MRFDIPname{5} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\TRUEFISP\NIST_241011_t600\bSSFP_DIP_Drop5_Epochs200'; %Nex600
MRFDIPname{6} = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\TRUEFISP\NIST_241011_t300\bSSFP_DIP_Drop5_Epochs200'; %Nex300

for j = 1: length(MRFDIPname)
    dataDIP{j} = load([MRFDIPname{j} '.mat']);
    t1map(:,:,j) = squeeze(dataDIP{j}.T1Iter(end,:,:));
    t2map(:,:,j) = squeeze(dataDIP{j}.T2Iter(end,:,:));
    m0map(:,:,j) = squeeze(dataDIP{j}.M0Iter(end,:,:));
    if isfield(dataDIP{j}, 'B0Iter')
        b0map(:,:,j) = squeeze(dataDIP{j}.B0Iter(end,:,:));
    else
        b0map(:,:,j) = zeros(size(squeeze(t1map(:,:,1))));
    end
    i=i+1;
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
for i = 1: length(MRFDIPname)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois);
    eval{i} = eval_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois, b0map(:,:,i));
end
%%
for i = 1:length(MRFDIPname)
    % [mean_est_t1, mean_est_t2] = plot_vs_NIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %
    [mean_est_t1(:,i), mean_est_t2(:,i)] = REGRNIST(NISTplot, t1map(:,:,i), t2map(:,:,i), rois); %

    % evaluation excluding short T1s and T2s
    th = 11;
    %regression
    [coeffs_t1(i,:), R2_t1(i)]=regression_fit_values(log(NISTplot(1:th, 1)), log(mean_est_t1(1:th,i)));
    [coeffs_t2(i,:), R2_t2(i)]=regression_fit_values(log(NISTplot(1:th, 2)), log(mean_est_t2(1:th,i)));
end
%% evaluation excluding short T1s and T2s


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

rel_error_tab_t1 = err_tab_t1./repmat(NISTplot(:,1),[1,6])*100;
rel_error_tab_t2 = err_tab_t2./repmat(NISTplot(:,2),[1,6])*100;

mean_rel_error_t1_th = mean(rel_error_tab_t1(1:th,:)); %excluding some vials
mean_rel_error_t2_th = mean(rel_error_tab_t2(1:th,:));
%% STD plots

f=figure('color','white');
f.Position = [300 300 1120 420];
b=bar(flipud(cov_t1)*100,'DisplayName','std_tab_t1');
b(1).FaceColor = [1 0 0]; b(2).FaceColor = [1 0.5 0]; b(3).FaceColor = [1 1 0];
b(4).FaceColor = [0 0 1]; b(5).FaceColor = [0 0.5 1]; b(6).FaceColor = [0 1 1];
xticklabels(flipud(NISTplot(:,1)))
title('COV T1 [%]')
legend('FISP-1000', 'FISP-600','FISP-300','trueFISP-1000','trueFISP-600','trueFISP-300', 'Location','northeast')

f=figure('color','white');
f.Position = [300 300 1120 420];
b=bar(flipud(cov_t2)*100,'DisplayName','std_tab_t2')
b(1).FaceColor = [1 0 0]; b(2).FaceColor = [1 0.5 0]; b(3).FaceColor = [1 1 0];
b(4).FaceColor = [0 0 1]; b(5).FaceColor = [0 0.5 1]; b(6).FaceColor = [0 1 1];
xticklabels(flipud(NISTplot(:,2)))
title('COV T2 [%]')
legend('FISP-1000', 'FISP-600','FISP-300','trueFISP-1000','trueFISP-600','trueFISP-300', 'Location','northeast')

%%
diff_FISP_t1cov(:,1) = cov_t1(:,1)-cov_t1(:,2);
diff_FISP_t1cov(:,2) = cov_t1(:,1)-cov_t1(:,3);
diff_FISP_t2cov(:,1) = cov_t2(:,1)-cov_t2(:,2);
diff_FISP_t2cov(:,2) = cov_t2(:,1)-cov_t2(:,3);

mFISP_t1 = mean(diff_FISP_t1cov,1);
mFISP_t2 = mean(diff_FISP_t2cov,1);


diff_trueFISP_t1cov(:,1) = cov_t1(:,4)-cov_t1(:,5);
diff_trueFISP_t1cov(:,2) = cov_t1(:,4)-cov_t1(:,6);
diff_trueFISP_t2cov(:,1) = cov_t2(:,4)-cov_t2(:,5);
diff_trueFISP_t2cov(:,2) = cov_t2(:,4)-cov_t2(:,6);

mtrueFISP_t1 = mean(diff_trueFISP_t1cov,1);
mtrueFISP_t2 = mean(diff_trueFISP_t2cov,1);
%%
% [h,p] = ttest(std_tab_t1(:,1) ,std_tab_t1(:,3))


%%


%%
s_vert = 100:300;
s_hor = 130:310;

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 3000];
climst2 = [0, 1500];
climstB0 = [-40 40];

for i = 1:length(MRFDIPname)
    f=figure('color','white');
    f.Position = [300 300 2400 420];
    ax1 = subplot(141); imagesc(t1map(s_vert,s_hor,i), climst1), colormap(ax1, T1colormap), colorbar, axis off
    ax2 = subplot(142); imagesc(t2map(s_vert,s_hor,i), climst2), colormap(ax2,T2colormap), colorbar, axis off
    mm = m0map(s_vert,s_hor,i);
    ax3 = subplot(143); imagesc(abs(mm)/max(abs(mm(:)))), colormap(ax3,gray), colorbar, axis off
    ax4 = subplot(144); imagesc(b0map(s_vert,s_hor,i),climstB0), colormap(ax4,gray), colorbar, axis off
end

%%
