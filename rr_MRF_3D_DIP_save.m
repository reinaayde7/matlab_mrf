% assembly of 3D DIP recon from 2D sequential DIP recon
clear all, close all, clc

parentFolder = uigetdir();
folderContents = dir(parentFolder);
folderNames = {folderContents([folderContents.isdir]).name};
folderNames = folderNames(~ismember(folderNames,{'.','..'}));
folderNames = natsort(folderNames);

% sTOT = 18; % original total number of partition
% j=6; % starting number of partition reconstructed
for i=1:length(folderNames)
    dataDIP{i} = load([parentFolder '\' folderNames{i} '\FISP_DIP_Drop5_Epochs400.mat']);
    % dataDIP{i} = load([parentFolder '\' folderNames{i} '\bSSFP_DIP_Drop5_Epochs400.mat']);
    
    t1map3d(:,:,i) = rot90(squeeze(dataDIP{i}.T1Iter(end-1,:,:)));
    t2map3d(:,:,i) = rot90(squeeze(dataDIP{i}.T2Iter(end-1,:,:)));
    m0map3d(:,:,i) = rot90(squeeze(dataDIP{i}.M0Iter(end-1,:,:)));
    if isfield(dataDIP{i}, 'B0Iter')
        b0map3d(:,:,i) = rot90(squeeze(dataDIP{i}.B0Iter(end,:,:)));
    else
        b0map3d(:,:,i) = zeros(size(t1map3d,1),size(t1map3d,2));
    end
    % j=j+1;
end
% 
% for i=j:sTOT
%     t1map3d(:,:,i) = zeros(size(t1map3d,1),size(t1map3d,2));
%     t2map3d(:,:,i) = zeros(size(t1map3d,1),size(t1map3d,2));
%     m0map3d(:,:,i) = zeros(size(t1map3d,1),size(t1map3d,2));
%     b0map3d(:,:,i) = zeros(size(t1map3d,1),size(t1map3d,2));
% end
%%
toSliceViewer = permute(dataDIP{1}.T1Iter(:,:,:),[2,3,1]);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar


%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 250]; %uplimit:1200 for NIST T2 at 0.55T 

toSliceViewer = t1map3d;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = t2map3d; 
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar
% 
toSliceViewer = abs(m0map3d);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

toSliceViewer = b0map3d;
figure,
sliceViewer(toSliceViewer, "DisplayRange", [-100, 100],'ScaleFactors', [1,1,1]), colorbar

%%
tempfilename = ['3D_MRF_DIP_' date '_DIP_3Dvolume'];
save([parentFolder '/' tempfilename], 't1map3d', 't2map3d', 'm0map3d', 'b0map3d')
