clear all, close all, clc
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500]; %1500/3000
climst2 = [0, 250]; %300/1200

folderdata = 'E:\scanner_data\twix_data\250221\';
nex600fa15 = load([folderdata '3D_MRF_fit_25-Feb-2025_meas_MID00547_brain_HV2_FISP_FA15_nex600_3Dfit_SVD.mat']);
nex1000fa15 = load([folderdata '3D_MRF_fit_25-Feb-2025_meas_MID00548_brain_HV2_FISP_FA15_nex1000_3Dfit_SVD.mat']);
nex600fa10 = load([folderdata '3D_MRF_fit_25-Feb-2025_meas_MID00549_brain_HV2_FISP_FA10_nex600_3Dfit_SVD.mat']);
nex1000fa10 = load([folderdata '3D_MRF_fit_25-Feb-2025_meas_MID00550_brain_HV2_FISP_FA10_nex1000_3Dfit_SVD.mat']);
%%
t1maps(:,:,:,1) = circshift(squeeze(nex600fa15.t1map3d),-2,3);
t1maps(:,:,:,2) = circshift(squeeze(nex1000fa15.t1map3d),-1,3);
t1maps(:,:,:,3) = squeeze(nex600fa10.t1map3d);
t1maps(:,:,:,4) = squeeze(nex1000fa10.t1map3d);

t2maps(:,:,:,1) = circshift(squeeze(nex600fa15.t2map3d),-2,3);
t2maps(:,:,:,2) = circshift(squeeze(nex1000fa15.t2map3d),-1,3);
t2maps(:,:,:,3) = squeeze(nex600fa10.t2map3d);
t2maps(:,:,:,4) = squeeze(nex1000fa10.t2map3d);

m0maps(:,:,:,1) = circshift(abs(squeeze(nex600fa15.m0map3d)),-2,3);
m0maps(:,:,:,2) = circshift(abs(squeeze(nex1000fa15.m0map3d)),-1,3);
m0maps(:,:,:,3) = abs(squeeze(nex600fa10.m0map3d));
m0maps(:,:,:,4) = abs(squeeze(nex1000fa10.m0map3d));

%% T1
d0 = squeeze(t1maps(:,:,:,1) - t1maps(:,:,:,1));
d1 = squeeze(t1maps(:,:,:,1) - t1maps(:,:,:,2));
d2 = squeeze(t1maps(:,:,:,1) - t1maps(:,:,:,3));
d3 = squeeze(t1maps(:,:,:,1) - t1maps(:,:,:,4));

line1 = cat(2,t1maps(:,:,:,1), t1maps(:,:,:,2), t1maps(:,:,:,3), t1maps(:,:,:,4));
line2 = cat(2,d0,d1,d2,d3);
toSliceViewer = cat(1,line1,line2);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%% T2
d0 = squeeze(t2maps(:,:,:,1) - t2maps(:,:,:,1));
d1 = squeeze(t2maps(:,:,:,1) - t2maps(:,:,:,2));
d2 = squeeze(t2maps(:,:,:,1) - t2maps(:,:,:,3));
d3 = squeeze(t2maps(:,:,:,1) - t2maps(:,:,:,4));

line1 = cat(2,t2maps(:,:,:,1), t2maps(:,:,:,2), t2maps(:,:,:,3), t2maps(:,:,:,4));
line2 = cat(2,d0,d1,d2,d3);
toSliceViewer = cat(1,line1,line2);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

%% M0
d0 = squeeze(m0maps(:,:,:,1) - m0maps(:,:,:,1));
d1 = squeeze(m0maps(:,:,:,1) - m0maps(:,:,:,2));
d2 = squeeze(m0maps(:,:,:,1) - m0maps(:,:,:,3));
d3 = squeeze(m0maps(:,:,:,1) - m0maps(:,:,:,4));

line1 = cat(2,m0maps(:,:,:,1), m0maps(:,:,:,2), m0maps(:,:,:,3), m0maps(:,:,:,4));
line2 = cat(2,d0,d1,d2,d3);
toSliceViewer = cat(1,line1,line2);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

%% import NIST data
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])

%% draw ROIs
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 200];
climstB0 = [-100 100];

codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))

% 
% slices = []
for i=1:4
    figure;
    imagesc(t1maps(:,:,33,i), climst1), colormap(T1colormap), colorbar
    axis image;
    rois{i} = get_rois;
end

%%
for r =1:length(rois)
    for m = 1:12
        [ix,iy] = find(rois{1,r}{1,m}==1);
        t1v{r,m}=reshape(squeeze(t1maps(ix,iy,33,r)),[size(ix,1)*size(iy,1),1]);
        
        t2v{r,m}=reshape(squeeze(t2maps(ix,iy,33,r)),[size(ix,1)*size(iy,1),1]);
    end
end

%% put rois together per tissue type

for meas = 1:4
    j=1;
    for r = 1:4:12
        v1{meas,j} = [t1v{meas,r}; t1v{meas,r+1}; t1v{meas,r+2}; t1v{meas,r+3}];
        t1meansM(meas,j) = mean(v1{meas,j});
        t1mediansM(meas,j) = median(v1{meas,j});
        t1stdsM(meas,j) = std(v1{meas,j});

        v2{meas,j} = [t2v{meas,r}; t2v{meas,r+1}; t2v{meas,r+2}; t2v{meas,r+3}];
        t2meansM(meas,j) = mean(v2{meas,j});
        t2mediansM(meas,j) = median(v2{meas,j});
        t2stdsM(meas,j) = std(v2{meas,j});

        j=j+1;

        
    end
    
end
% gm_m1 = [t1v{1,1}; t1v{1,2}; t1v{1,3}; t1v{1,4}];
% gm_m2 = [t1v{2,1}; t1v{2,2}; t1v{2,3}; t1v{2,4}];
% gm_m3 = [t1v{3,1}; t1v{3,2}; t1v{3,3}; t1v{3,4}];
% gm_m4 = [t1v{4,1}; t1v{4,2}; t1v{4,3}; t1v{4,4}];
