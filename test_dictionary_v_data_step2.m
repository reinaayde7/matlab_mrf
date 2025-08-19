clear all, close all, clc

%%
disp('open processed data')
RawData.folder = 'E:\scanner_data\twix_data\241011\';
RawData.nameTRUEFISP = 'DataStored_img_rawdatainfo_21-Oct-2024_meas_MID01027_trueFISPv2.mat';
RawData.nameFISP = 'DataStored_img_rawdatainfo_21-Oct-2024_meas_MID01025_FISP.mat';
RawData.b0map = 'meas_MID01030_FID112465_b0_mapping_rudy_OFFRES';


dFISP = load([RawData.folder RawData.nameFISP]);
dTRUEFISP = load([RawData.folder RawData.nameTRUEFISP]);
b0map = load([RawData.folder RawData.b0map]);

%% draw ROIs: NIST 1-14 ordering
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(dFISP.proxy_img));
axis image;
rois = get_rois;

%% average signal in each ROI
NumFrames = size(dFISP.image_combined,3);
avg_fp = zeros(size(rois,2),NumFrames);
for roiIDX = 1: size(rois,2)
    [ix,iy] = find(rois{roiIDX}==1);
    fp= reshape(dFISP.image_combined(ix,iy,:),[size(ix,1)*size(iy,1),NumFrames]);
    dFISP.avg_fp(roiIDX,:) = mean(fp);
    tfp= reshape(dTRUEFISP.image_combined(ix,iy,:),[size(ix,1)*size(iy,1),NumFrames]);
    dTRUEFISP.avg_fp(roiIDX,:) = mean(tfp);
end

%% average b0 values in each vial
for roiIDX = 1: size(rois,2)
    [ix,iy] = find(rois{roiIDX}==1);
    b0values= reshape(b0map.LPF(ix,iy,:),[size(ix,1)*size(iy,1),1]);
    b0map.avg_b0(roiIDX) = mean(b0values);
end


%% save for comparison
tempfilename = ['SingleVialAnalysis_trueFISPv2vFISP_withB0map_' date];
RawData.folder = 'E:\scanner_data\twix_data\241011\';
save([RawData.folder tempfilename], 'dFISP', 'dTRUEFISP', 'b0map')