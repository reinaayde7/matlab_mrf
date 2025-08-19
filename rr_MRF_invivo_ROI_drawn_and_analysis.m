% rr MRF 3D plot and gifs
clear all, close all, clc

%%

MRFdir = ['E:\POSTDOC_UoM\08_Project_MRF\ISMRM_25_datasets\'];

MRFname = { '3D_MRF_DIP_02-Nov-2024_brain_FISP_v1'; 
'3D_MRF_DIP_02-Nov-2024_brain_FISP_v2';
'3D_MRF_DIP_02-Nov-2024_brain_FISP_v3';
'3D_MRF_DIP_02-Nov-2024_brain_FISP_v4';
'3D_MRF_fit_23-Oct-2024_meas_MID00341_brain_FISP_v1';
'3D_MRF_fit_01-Nov-2024_meas_MID00051_brain_FISP_v2';
'3D_MRF_fit_01-Nov-2024_meas_MID00274_brain_FISP_v3';
'3D_MRF_fit_01-Nov-2024_meas_MID00287_brain_FISP_v4';
'3D_MRF_DIP_02-Nov-2024_brain_TRUEFISP_v1';
'3D_MRF_DIP_02-Nov-2024_brain_TRUEFISP_v2';
'3D_MRF_DIP_02-Nov-2024_brain_TRUEFISP_v3';
'3D_MRF_DIP_02-Nov-2024_brain_TRUEFISP_v4';
};
for i = 1: length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,:,i) = squeeze(data{i}.t1map3d);
    t2map(:,:,:,i) = squeeze(data{i}.t2map3d);
    m0map(:,:,:,i) = squeeze(data{i}.m0map3d);
    b0map(:,:,:,i) = rot90(squeeze(data{i}.b0map3d));
end


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
for i=1:4
    figure;
    imagesc(t1map(:,:,slice,i), climst1), colormap(T1colormap), colorbar
    axis image;
    rois{i} = get_rois;
end
%%
for j = 1:4:length(MRFname)
    for r =1:length(rois)
        for m = 1:9
            [ix,iy] = find(rois{1,r}{1,m}==1);
            fp=reshape(squeeze(t1map(ix,iy,slice,j+r-1)),[size(ix,1)*size(iy,1),1]);
            avgV(r+j-1,m) = mean(fp);
            stdV(r+j-1,m) = std(fp);
        end
    end
end



%% print plot

idx1 = 130:290;
idy1 = 55:260;

idx2 = 120:290;
idy2 = 40:260;

idx3 = 125:270;
idy3 = 85:270;

idx4 = 115:280;
idy4 = 85:270;

slice = 11;
fig=figure; 
fig.Position = [0 0 1100 800];
t = tiledlayout(3, 4, "TileSpacing", "none");

%-----------------T1 maps------------
colormapP = T1colormap;
climst = climst1;

%v1
ax = nexttile;
ImageAux = t1map(idy1,idx1,slice-2,9);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
text(0.5, 1.1, 'HV #1', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
%v2
ax = nexttile;
ImageAux = t1map(idy2,idx2,slice,10);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(0.5, 1.1, 'HV #2', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
%v3
ax = nexttile;
ImageAux = t1map(idy3,idx3,slice+1,11);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(0.5, 1.1, 'HV #3', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
%v4
ax = nexttile;
ImageAux = t1map(idy4,idx4,slice+1,12);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(0.5, 1.1, 'HV #4', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;

%-----------------T2 maps------------
colormapP = T2colormap;
climst = climst2;

%v1
ax = nexttile;
ImageAux = t2map(idy1,idx1,slice-2,9);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
%v2
ax = nexttile;
ImageAux = t2map(idy2,idx2,slice,10);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
%v3
ax = nexttile;
ImageAux = t2map(idy3,idx3,slice+1,11);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
%v4
ax = nexttile;
ImageAux = t2map(idy4,idx4,slice+1,12);
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;

%-----------------M0 maps------------
colormapP = 'gray';
climst = [0, 1];

%v1
ax = nexttile;
mm0 = squeeze(abs(m0map(idy1,idx1,slice-2,9)));
m0plot = mm0./max(mm0(:));
ImageAux = m0plot;
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);

%v2
ax = nexttile;
mm0 = squeeze(abs(m0map(idy2,idx2,slice,10)));
m0plot = mm0./max(mm0(:));
ImageAux = m0plot;
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
%v3
ax = nexttile;
mm0 = squeeze(abs(m0map(idy3,idx3,slice+1,11)));
m0plot = mm0./max(mm0(:));
ImageAux = m0plot;
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
%v4
ax = nexttile;
mm0 = squeeze(abs(m0map(idy4,idx4,slice+1,12)));
m0plot = mm0./max(mm0(:));
ImageAux = m0plot;
imagesc(ImageAux, climst), colormap(ax, colormapP)
% axis equal
pbaspect([1 1.1 1])
axis off
cb0 = colorbar(ax);
cb0.Position(4) = cb0.Position(4)*0.9;