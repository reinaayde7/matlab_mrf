% rr MRF 3D plot and gifs
clear all, close all, clc

%%

% HV5
MRFdir = ['E:\scanner_data\twix_data\250305\'];
MRFname = { '3D_MRF_fit_06-Mar-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD',
    '3D_MRF_fit_14-Mar-2025_meas_MID00017_brain_HV5_bSSFP_FA15_nex600_w3000_3Dfit_SVD',
    '3D_MRF_DIP_07-Mar-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD',
    '3D_MRF_DIP_13-Mar-2025_meas_MID00017_brain_HV5_bSSFP_FA15_nex600_w3000_3Dfit_SVD'
    };

% HV4
% MRFdir = ['E:\scanner_data\twix_data\250228\'];
% MRFname = { '3D_MRF_fit_05-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD',
%     '3D_MRF_fit_19-Mar-2025_meas_MID00026_brain_HV4_bSSFP_FA15_nex600_w3000_3Dfit_SVD',
%     '3D_MRF_DIP_07-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD',
%     '3D_MRF_DIP_13-Mar-2025_meas_MID00026_brain_HV4_bSSFP_FA15_nex600_3Dfit_SVD'
%     };

saveFolder = 'plots_LRvDIP_HV5_FISPvpSSFP_FA15_nex600_withabsB0\';
mkdir([MRFdir saveFolder])
save_filenames = 'plots_LRvDIP_HV5_FISPvpSSFP_FA15_nex600';
%%
for i = 1: length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    t1map(:,:,:,i) = squeeze(data{i}.t1map3d);
    t2map(:,:,:,i) = squeeze(data{i}.t2map3d);
    m0map(:,:,:,i) = squeeze(data{i}.m0map3d);
    b0map(:,:,:,i) = abs(squeeze(data{i}.b0map3d));
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
idx1 = 5:235; %vertical 
idx2 = 45:275; %horizontal
%v2 zoom in
idx1z = 30:100; %vertical
idx2z = 150:220; %horizontal

toSliceViewer = squeeze(t1map(idx1z,idx2z,:,1));

figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%%
toplot = cat(2,squeeze(t1map(idx1,idx2,:,1)),squeeze(t1map(idx1,idx2,:,2)),squeeze(t1map(idx1,idx2,:,3)),squeeze(t1map(idx1,idx2,:,4)));
figure,
sliceViewer(toplot, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar



%%
toplot = cat(2,squeeze(t2map(idx1,idx2,:,1)),squeeze(t2map(idx1,idx2,:,2)),squeeze(t2map(idx1,idx2,:,3)),squeeze(t2map(idx1,idx2,:,4)));
figure,
sliceViewer(toplot, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar


%%

%single slice scaling for intensity
for s=1:48
    c1 = abs(squeeze(m0map(:,:,s,1)));
    c1m0(:,:,s) = c1/max(c1(:))*1.8;
    c2 = abs(squeeze(m0map(:,:,s,2)));
    c2m0(:,:,s) = c2/max(c2(:))*1.4;
    c3 = abs(squeeze(m0map(:,:,s,3)));
    c3m0(:,:,s) = c3/max(c3(:));
    c4 = abs(squeeze(m0map(:,:,s,4)));
    c4m0(:,:,s) = c4/max(c4(:));
end

% %global scaling for intensity
% c1m0 = abs(squeeze(m0map(:,:,:,1)));
% c1m0 = c1m0/max(c1m0(:))*1.8;
% c2m0 = abs(squeeze(m0map(:,:,:,2)));
% c2m0 = c2m0/max(c2m0(:))*1.8;
% c3m0 = abs(squeeze(m0map(:,:,:,3)));
% c3m0 = c3m0/max(c3m0(:));
% c4m0 = abs(squeeze(m0map(:,:,:,4)));
% c4m0 = c4m0/max(c4m0(:))*4;

mm0(:,:,:,1) = c1m0;
mm0(:,:,:,2) = c2m0;
mm0(:,:,:,3) = c1m0-c2m0;

mm0(:,:,:,1) = c1m0;
mm0(:,:,:,2) = c2m0;
mm0(:,:,:,3) = c3m0;
mm0(:,:,:,4) = c4m0;


climM0 = [0,0.7]; %0.6
toplot = cat(2,c1m0(idx1,idx2,:),c2m0(idx1,idx2,:),c3m0(idx1,idx2,:),c4m0(idx1,idx2,:));
figure,
sliceViewer(toplot, "DisplayRange", climM0, 'ScaleFactors', [1,1,1]), colorbar

%% mask calculation
figure, histogram(c3m0)
th = 0.04;
for s = 1:size(c2m0,3)
    M0s = c3m0(:,:,s);
    cluster_masks(:,:,s) = M0s>th;
end

figure,
sliceViewer(cluster_masks), colorbar

%% masking
for s = 1:size(c2m0,3)
    for j =1:4
        t1map(:,:,s,j) = t1map(:,:,s,j).*cluster_masks(:,:,s); 
        t2map(:,:,s,j) = t2map(:,:,s,j).*cluster_masks(:,:,s); 
        mm0(:,:,s,j) = mm0(:,:,s,j).*cluster_masks(:,:,s); 
    end
end


%% 1.2 video generator
% design from ISMRM25 where 3 different reconstructions where compared

filename =[MRFdir saveFolder save_filenames];
indexes = [10:34];
for n = 1:numel(indexes)
    slice = indexes(n);


    Image = cat(4,squeeze(t1map(idx1,idx2,slice,:)),squeeze(t2map(idx1,idx2,slice,:)),squeeze(mm0(idx1,idx2,slice,:)));
    ImageZ = cat(4,squeeze(t1map(idx1z,idx2z,slice,:)),squeeze(t2map(idx1z,idx2z,slice,:)),squeeze(mm0(idx1z,idx2z,slice,:)));
    
    fig=figure; 
    fig.Position = [0 0 1700 900];
    fig.Color = 'white';
    t = tiledlayout(4, 8, "TileSpacing", "none");
    
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
    text(0.5, 1.1, 'LR-FISP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    rectangle('Position', [105, 25, 70, 70], 'EdgeColor', 'y', 'LineWidth', 2);
    ax = nexttile;
    ImageAux = Image(:,:,2,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'LR-pSSFP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = Image(:,:,3,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'DIP-FISP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = Image(:,:,4,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'DIP-pSSFP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    

    ax = nexttile;
    ImageAux = ImageZ(:,:,1,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'LR-FISP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,2,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'LR-pSSFP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,3,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'DIP-FISP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,4,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'DIP-pSSFP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
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
    ImageAux = Image(:,:,4,2);
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
     ax = nexttile;
    ImageAux = ImageZ(:,:,4,2);
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
    ImageAux = Image(:,:,4,3);
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
    ax = nexttile;
    ImageAux = ImageZ(:,:,4,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    cb2.FontSize = 12;
    % cb2.Position(1) = cb2.Position(1)+0.03;

     
    ax = nexttile;axis off;
    text(-0.1, 0.5, '|B0| [Hz]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    ax = nexttile;
    bb = rot90(b0map(:,:,slice,2),1);
    ImageAux = bb(idx1,idx2).*cluster_masks(idx1,idx2,slice); 
    imagesc(ImageAux, [0, 100]), colormap(ax, gray)
    axis off
 
    ax = nexttile;  
    axis off;  
    ax = nexttile;
    ImageAux = b0map(idx1,idx2,slice,4).*cluster_masks(idx1,idx2,slice); 
    imagesc(ImageAux, [0, 100]), colormap(ax, gray)
    axis off

    ax = nexttile;axis off;
    ax = nexttile;
    bb = rot90(b0map(:,:,slice,2),1);
    ImageAux = bb(idx1z,idx2z).*cluster_masks(idx1z,idx2z,slice); 
    imagesc(ImageAux, [0, 100]), colormap(ax, gray)
    axis off

    ax = nexttile;axis off;
    ax = nexttile;
    ImageAux = b0map(idx1z,idx2z,slice,4).*cluster_masks(idx1z,idx2z,slice); 
    imagesc(ImageAux, [0, 100]), colormap(ax, gray)
    axis off
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    cb2.FontSize = 12;
    
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

