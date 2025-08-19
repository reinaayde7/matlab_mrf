% rr MRF 3D plot and gifs
clear all, close all, clc

%%

%volunteer 241018 3D LR recon 
% MRFdir = ['E:\POSTDOC_UoM\08_Project_MRF\ISMRM_25_datasets\'];
MRFdir = ['E:\scanner_data\twix_data\250203\'];
MRFname = { '3D_MRF_fit_07-Apr-2025_meas_MID00188_brain_HV1_FISP_FA15_nex600_w3000_3Dfit_SVD_LR_SLLR',
    '3D_MRF_DIP_12-Feb-2025_MID00188_DIP_3Dvolume'
    };

saveFolder = 'plots_LRvSLLRvDIP_HV1_FISP_FA15_nex600\';
mkdir([MRFdir saveFolder])
save_filenames = 'plots_LRvSLLRvDIP_HV1_FISP_FA15_nex600';
%%
k=1;
for i = 1:1%length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,:,k) = rot90(squeeze(data{i}.t1map3d));
    t2map(:,:,:,k) = rot90(squeeze(data{i}.t2map3d));
    m0map(:,:,:,k) = rot90(squeeze(data{i}.m0map3d));
    b0map(:,:,:,k) = rot90(squeeze(data{i}.b0map3d));
 
    k=k+1;
    t1map(:,:,:,k) = rot90(squeeze(data{i}.t1map3dLLR));
    t2map(:,:,:,k) = rot90(squeeze(data{i}.t2map3dLLR));
    m0map(:,:,:,k) = rot90(squeeze(data{i}.m0map3dLLR));
    b0map(:,:,:,k) = rot90(squeeze(data{i}.b0map3dLLR));

    k=k+1;
end

for i = 2:2%length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,:,k) = squeeze(data{i}.t1map3d);
    t2map(:,:,:,k) = squeeze(data{i}.t2map3d);
    m0map(:,:,:,k) = squeeze(data{i}.m0map3d);
    b0map(:,:,:,k) = squeeze(data{i}.b0map3d);
 

    k=k+1;
end

%% import acquired b0 map
% b0name = 'meas_MID00343_FID113798_b0_mapping_rudy_200x200_st3_5_OFFRES';
% load([MRFdir b0name '.mat']);
% LPF = rot90(LPF);
%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 250]; %uplimit:1200 for NIST T2 at 0.55T 

%% 1.1 test display (make sure they are a square to keep proportion and res)
%v2 full
idx1 = 30:235; %vertical 
idx2 = 35:240; %horizontal
%v2 zoom in
idx1z = 50:120; %vertical
idx2z = 130:200; %horizontal

toSliceViewer = squeeze(t1map(idx1z,idx2z,:,1));

figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%%
% t1map(:,:,:,3) = squeeze(t1map(:,:,:,1))-squeeze(t1map(:,:,:,2));
toplot = cat(2,squeeze(t1map(idx1,idx2,:,1)),squeeze(t1map(idx1,idx2,:,2)),squeeze(t1map(idx1,idx2,:,3)));
figure,
sliceViewer(toplot, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar



%%
% t2map(:,:,:,3) = squeeze(t2map(:,:,:,1))-squeeze(t2map(:,:,:,2));
toplot = cat(2,squeeze(t2map(idx1,idx2,:,1)),squeeze(t2map(idx1,idx2,:,2)),squeeze(t2map(idx1,idx2,:,3)));
figure,
sliceViewer(toplot, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar


%%
c1m0 = abs(squeeze(m0map(:,:,:,1)));
c1m0 = c1m0/max(c1m0(:))*1.8;
c2m0 = abs(squeeze(m0map(:,:,:,2)));
c2m0 = c2m0/max(c2m0(:))*1.8;
c3m0 = abs(squeeze(m0map(:,:,:,3)));
c3m0 = c3m0/max(c3m0(:));

mm0(:,:,:,1) = c1m0;
mm0(:,:,:,2) = c2m0;
% mm0(:,:,:,3) = c1m0-c2m0;
mm0(:,:,:,3) = c3m0;

climM0 = [0,0.6]; %0.6
% toplot = cat(2,c1m0(idx1,idx2,:),c2m0(idx1,idx2,:),c1m0(idx1,idx2,:)-c2m0(idx1,idx2,:));
toplot = cat(2,c1m0(idx1,idx2,:),c2m0(idx1,idx2,:),c3m0(idx1,idx2,:));
figure,
sliceViewer(toplot, "DisplayRange", climM0, 'ScaleFactors', [1,1,1]), colorbar

%% mask calculation
figure, histogram(c3m0)
th = 0.03;
for s = 1:size(c3m0,3)
    M0s = c3m0(:,:,s);
    cluster_masks(:,:,s) = M0s>th;
end

figure,
sliceViewer(cluster_masks), colorbar

%% masking
for s = 1:size(c3m0,3)
    for j =1:3
        t1map(:,:,s,j) = t1map(:,:,s,j).*cluster_masks(:,:,s); 
        t2map(:,:,s,j) = t2map(:,:,s,j).*cluster_masks(:,:,s); 
        mm0(:,:,s,j) = mm0(:,:,s,j).*cluster_masks(:,:,s); 
    end
end


%% 1.2 video generator
% design from ISMRM25 where 3 different reconstructions where compared

filename =[MRFdir saveFolder save_filenames];
indexes = [2:33];
for n = 1:numel(indexes)
    slice = indexes(n);


    Image = cat(4,squeeze(t1map(idx1,idx2,slice,:)),squeeze(t2map(idx1,idx2,slice,:)),squeeze(mm0(idx1,idx2,slice,:)));
    ImageZ = cat(4,squeeze(t1map(idx1z,idx2z,slice,:)),squeeze(t2map(idx1z,idx2z,slice,:)),squeeze(mm0(idx1z,idx2z,slice,:)));
    
    fig=figure; 
    fig.Position = [0 0 1600 800];
    fig.Color = 'white';
    t = tiledlayout(3, 6, "TileSpacing", "none");
    
    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
    
    ax = nexttile;
    
    ImageAux = Image(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % axis equal
    % pbaspect([1 length(idx1)/length(idx2) 1])
    axis off
    text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    text(0.5, 1.1, 'Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    rectangle('Position', [95, 20, 70, 70], 'EdgeColor', 'y', 'LineWidth', 2);
    ax = nexttile;
    ImageAux = Image(:,:,2,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'Spatial-Local-Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = Image(:,:,3,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'Deep-Image-Prior', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'Spatial-Local-Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'Deep-Image-Prior', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    cb0 = colorbar(ax);
    % cb0.Position(1) = cb0.Position(1)+0.05;
    cb0.Position(4) = cb0.Position(4)*0.9;
    cb0.FontSize = 12;
    
    
    %T2 maps
    colormapP = T2colormap;
    climst = climst2;
    
    ax = nexttile;
    ImageAux = Image(:,:,1,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    ax = nexttile;
    ImageAux = Image(:,:,2,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    ax = nexttile;
    ImageAux = Image(:,:,3,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    cb1.FontSize = 12;
    % cb1.Position(1) = cb1.Position(1)+0.037;
    
    
    %M0 maps
    colormapP = 'gray';
    climst = climM0;
    ax = nexttile;
    ImageAux = Image(:,:,1,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    ax = nexttile;
    ImageAux = Image(:,:,2,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    ax = nexttile;
    ImageAux = Image(:,:,3,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    cb2.FontSize = 12;
    % cb2.Position(1) = cb2.Position(1)+0.03;
    
    exportgraphics(gcf,[MRFdir saveFolder save_filenames '_s_' num2str(slice) '.tif'], 'Resolution',300)
    %gif saves
    timer=1;
    frame = getframe(fig);
    im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
            imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',timer);
        else
            imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',timer);
        end
        
end

%% 1.2 differences
% design from ISMRM25 where 3 different reconstructions where compared

filename =[MRFdir saveFolder 'DIFF_' save_filenames];
indexes = [2:33];
for n = 1:numel(indexes)
    slice = indexes(n);


    Image = cat(4,squeeze(t1map(idx1,idx2,slice,:)),squeeze(t2map(idx1,idx2,slice,:)),squeeze(mm0(idx1,idx2,slice,:)));
    ImageZ = cat(4,squeeze(t1map(idx1z,idx2z,slice,:)),squeeze(t2map(idx1z,idx2z,slice,:)),squeeze(mm0(idx1z,idx2z,slice,:)));
    
    fig=figure; 
    fig.Position = [0 0 870 800];
    fig.Color = 'white';
    t = tiledlayout(3, 3, "TileSpacing", "none");
    
    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
   
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,1)-ImageZ(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    text(0.5, 1.1, 'Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,1)-ImageZ(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'Spatial-Local-Low-Rank', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,1)-ImageZ(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'Deep-Image-Prior', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    cb0 = colorbar(ax);
    % cb0.Position(1) = cb0.Position(1)+0.05;
    cb0.Position(4) = cb0.Position(4)*0.9;
    cb0.FontSize = 12;
    
    
    %T2 maps
    colormapP = T2colormap;
    climst = climst2;
    
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,2)-ImageZ(:,:,1,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,2)-ImageZ(:,:,1,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,2)-ImageZ(:,:,1,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    cb1.FontSize = 12;
    % cb1.Position(1) = cb1.Position(1)+0.037;
    
    
    %M0 maps
    colormapP = 'gray';
    climst = climM0;
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1,3)-ImageZ(:,:,1,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,3)-ImageZ(:,:,1,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,3)-ImageZ(:,:,1,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    cb2.FontSize = 12;
    % cb2.Position(1) = cb2.Position(1)+0.03;
    
    exportgraphics(gcf,[MRFdir saveFolder 'DIFF_' save_filenames '_s_' num2str(slice) '.tif'], 'Resolution',300)
    %gif saves
    timer=1;
    frame = getframe(fig);
    im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
            imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',timer);
        else
            imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',timer);
        end
        
end