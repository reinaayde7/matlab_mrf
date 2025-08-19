% rr MRF analysis
%considering one slice only!
clear all, close all, clc

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);

saveFolder = 'plots_denoisingMULTIPLES_HV3_FISP_FA15_nex600\';
mkdir([MRFdir saveFolder])
save_filenames = 'plots_denoisingMULTIPLES_FISP_FA15_nex600';
%%
MRFdir = ['E:\scanner_data\twix_data\250221\'];
MRFname = { 

    '3D_MRF_fit_27-Feb-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_3Dfit_SVD';
    '3D_MRF_fit_27-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_bm3d';
    '3D_MRF_fit_27-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_bm4d';
    '3D_MRF_fit_27-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_MPPCA';

    };

for i = 1: length(MRFname)
    data{i} = load([MRFdir MRFname{i} '.mat']);
    if i < 2
        t1map(:,:,:,i) = data{i}.t1map3d;
        t2map(:,:,:,i) = data{i}.t2map3d;
        m0map(:,:,:,i) = data{i}.m0map3d;
    else
        t1map(:,:,:,i) = rot90(data{i}.t1map3d);
        t2map(:,:,:,i) = rot90(data{i}.t2map3d);
        m0map(:,:,:,i) = rot90(data{i}.m0map3d);
    end
end


%% dip import

MRFname = { 
    '3D_MRF_DIP_03-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_DIP_3Dvolume_HV3';
    };

 data{i+1} = load([MRFdir MRFname{1} '.mat']);

 t1map(:,:,:,i+1) = data{i+1}.t1map3d;
 t2map(:,:,:,i+1) = data{i+1}.t2map3d;
 m0map(:,:,:,i+1) = data{i+1}.m0map3d;

%%

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 1500];
climst2 = [0, 250];

%v2 full
idx1 = 10:225; %vertical 
idx2 = 50:265; %horizontal
%v2 zoom in
idx1z = 40:110; %vertical
idx2z = 150:220; %horizontal

t1map(:,:,:,5) = circshift(t1map(:,:,:,5),-1,1);
toSliceViewer = squeeze(t1map(idx1,idx2,:,1));


figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%%
c1m0 = abs(squeeze(m0map(:,:,:,1)));
c1m0 = c1m0/max(c1m0(:))*1.8;
c2m0 = abs(squeeze(m0map(:,:,:,2)));
c2m0 = c2m0/max(c2m0(:))*1.8;
c3m0 = abs(squeeze(m0map(:,:,:,3)));
c3m0 = c3m0/max(c3m0(:))*1.8;
c4m0 = abs(squeeze(m0map(:,:,:,4)));
c4m0 = c4m0/max(c4m0(:))*1.8;
c5m0 = abs(squeeze(m0map(:,:,:,5)));
c5m0 = c5m0/max(c5m0(:));

mm0(:,:,:,1) = c1m0;
mm0(:,:,:,2) = c2m0;
mm0(:,:,:,3) = c3m0;
mm0(:,:,:,4) = c4m0;
mm0(:,:,:,5) = c5m0;

climM0 = [0,0.6]; %0.6
toplot = cat(3,c1m0(idx1,idx2,:),c2m0(idx1,idx2,:),c3m0(idx1,idx2,:),c4m0(idx1,idx2,:),c5m0(idx1,idx2,:));
figure,
sliceViewer(toplot, "DisplayRange", climM0, 'ScaleFactors', [1,1,1]), colorbar

%% mask calculation
figure, histogram(c5m0)
th = 0.04;
for s = 1:size(c5m0,3)
    M0s = c5m0(:,:,s);
    cluster_masks(:,:,s) = M0s>th;
end

figure,
sliceViewer(cluster_masks), colorbar


%% masking
for s = 1:size(c2m0,3)
    for j =1:5
        t1map(:,:,s,j) = t1map(:,:,s,j).*cluster_masks(:,:,s); 
        t2map(:,:,s,j) = t2map(:,:,s,j).*cluster_masks(:,:,s); 
        mm0(:,:,s,j) = mm0(:,:,s,j).*cluster_masks(:,:,s); 
    end
end



%% overall plotting

filename =[MRFdir saveFolder save_filenames];
indexes = [10:38];
for n = 1:numel(indexes)
    slice = indexes(n);

    Image = cat(3,squeeze(t1map(idx1,idx2,slice,1)),squeeze(t1map(idx1,idx2,slice,2)),squeeze(t1map(idx1,idx2,slice,3)),squeeze(t1map(idx1,idx2,slice,4)),squeeze(t1map(idx1,idx2,slice,5)));
    ImageZ = cat(3,squeeze(t1map(idx1z,idx2z,slice,1)),squeeze(t1map(idx1z,idx2z,slice,2)),squeeze(t1map(idx1z,idx2z,slice,3)),squeeze(t1map(idx1z,idx2z,slice,4)),squeeze(t1map(idx1z,idx2z,slice,5)));
    
    fig=figure; 
    fig.Position = [0 0 1800 500];
    fig.Color = 'white';
    t = tiledlayout(3, 10, "TileSpacing", "none");
    
    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % axis equal
    % pbaspect([1 length(idx1)/length(idx2) 1])
    axis off
    text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    text(0.5, 1.1, 'LR', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    rectangle('Position', [100, 30, 70, 70], 'EdgeColor', 'y', 'LineWidth', 2);
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'BM3D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'BM4D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'MPPCA', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'DIP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'LR', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'BM3D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'BM4D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'MPPCA', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'DIP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    cb0 = colorbar(ax);
    % cb0.Position(1) = cb0.Position(1)+0.05;
    cb0.Position(4) = cb0.Position(4)*0.9;
    cb0.FontSize = 12;
    
    
    %T2 maps
    Image = cat(3,squeeze(t2map(idx1,idx2,slice,1)),squeeze(t2map(idx1,idx2,slice,2)),squeeze(t2map(idx1,idx2,slice,3)),squeeze(t2map(idx1,idx2,slice,4)),squeeze(t2map(idx1,idx2,slice,5)));
    ImageZ = cat(3,squeeze(t2map(idx1z,idx2z,slice,1)),squeeze(t2map(idx1z,idx2z,slice,2)),squeeze(t2map(idx1z,idx2z,slice,3)),squeeze(t2map(idx1z,idx2z,slice,4)),squeeze(t2map(idx1z,idx2z,slice,5)));
    
    colormapP = T2colormap;
    climst = climst2;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    cb1.FontSize = 12;
    % cb1.Position(1) = cb1.Position(1)+0.037;
    
    
    %M0 maps
    Image = cat(3,squeeze(mm0(idx1,idx2,slice,1)),squeeze(mm0(idx1,idx2,slice,2)),squeeze(mm0(idx1,idx2,slice,3)),squeeze(mm0(idx1,idx2,slice,4)),squeeze(mm0(idx1,idx2,slice,5)));
    ImageZ = cat(3,squeeze(mm0(idx1z,idx2z,slice,1)),squeeze(mm0(idx1z,idx2z,slice,2)),squeeze(mm0(idx1z,idx2z,slice,3)),squeeze(mm0(idx1z,idx2z,slice,4)),squeeze(mm0(idx1z,idx2z,slice,5)));
    colormapP = 'gray';
    climst = climM0;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
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

%% differences

filename =[MRFdir saveFolder save_filenames];
indexes = [10:38];
for n = 1:numel(indexes)
    slice = indexes(n);

    for dentyp = 1:5
        d(:,:,dentyp) = t1map(:,:,slice,1)- t1map(:,:,slice,dentyp);
    end
    
    ii = t1map(:,:,slice,:)>2000; %filter out CSF
    d(ii) = 0;%+10*rand(1);
    
    Image = cat(3,squeeze(d(idx1,idx2,1)),squeeze(d(idx1,idx2,2)),squeeze(d(idx1,idx2,3)),squeeze(d(idx1,idx2,4)),squeeze(d(idx1,idx2,5)));
    ImageZ = cat(3,squeeze(d(idx1z,idx2z,1)),squeeze(d(idx1z,idx2z,2)),squeeze(d(idx1z,idx2z,3)),squeeze(d(idx1z,idx2z,4)),squeeze(d(idx1z,idx2z,5)));
    
    fig=figure; 
    fig.Position = [0 0 1800 500];
    fig.Color = 'white';
    t = tiledlayout(3, 10, "TileSpacing", "none");
    
    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % axis equal
    % pbaspect([1 length(idx1)/length(idx2) 1])
    axis off
    text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    text(0.5, 1.1, 'LR', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    rectangle('Position', [105, 30, 70, 70], 'EdgeColor', 'y', 'LineWidth', 2);
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'BM3D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'BM4D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'MPPCA', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    text(0.5, 1.1, 'DIP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'LR', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'BM3D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'BM4D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'MPPCA', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    text(0.5, 1.1, 'DIP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    
    cb0 = colorbar(ax);
    % cb0.Position(1) = cb0.Position(1)+0.05;
    cb0.Position(4) = cb0.Position(4)*0.9;
    cb0.FontSize = 12;
    
    %T2 maps
    for dentyp = 1:5
        d(:,:,dentyp) = t2map(:,:,slice,1)-t2map(:,:,slice,dentyp);
    end
    
    % ii = t1map>1000; %filter out CSF
    d(ii) = 0;%+20*rand(1);
    
    Image = cat(3,squeeze(d(idx1,idx2,1)),squeeze(d(idx1,idx2,2)),squeeze(d(idx1,idx2,3)),squeeze(d(idx1,idx2,4)),squeeze(d(idx1,idx2,5)));
    ImageZ = cat(3,squeeze(d(idx1z,idx2z,1)),squeeze(d(idx1z,idx2z,2)),squeeze(d(idx1z,idx2z,3)),squeeze(d(idx1z,idx2z,4)),squeeze(d(idx1z,idx2z,5)));
    
    colormapP = T2colormap;
    climst = climst2;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    cb1.FontSize = 12;
    % cb1.Position(1) = cb1.Position(1)+0.037;
    
    
    %M0 maps
    for dentyp = 1:5
        d(:,:,dentyp) = mm0(:,:,slice,1)-mm0(:,:,slice,dentyp);
    end
    d(ii) = 0;
    
    Image = cat(3,squeeze(d(idx1,idx2,1)),squeeze(d(idx1,idx2,2)),squeeze(d(idx1,idx2,3)),squeeze(d(idx1,idx2,4)),squeeze(d(idx1,idx2,5)));
    ImageZ = cat(3,squeeze(d(idx1z,idx2z,1)),squeeze(d(idx1z,idx2z,2)),squeeze(d(idx1z,idx2z,3)),squeeze(d(idx1z,idx2z,4)),squeeze(d(idx1z,idx2z,5)));
    colormapP = 'gray';
    climst = climM0;
    
    ax = nexttile;
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = Image(:,:,5);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % pbaspect([1 1.2 1])
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,4);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    axis off
    
    ax = nexttile;
    ImageAux = ImageZ(:,:,5);
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
