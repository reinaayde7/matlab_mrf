% rr MRF 3D plot and gifs
clear all, close all, clc

%%

%volunteer 241018 3D LR recon 
% MRFdir = ['E:\POSTDOC_UoM\08_Project_MRF\ISMRM_25_datasets\'];
MRFdir = {'E:\scanner_data\twix_data\250305\',
    'E:\scanner_data\twix_data\250228\',
    'E:\scanner_data\twix_data\250221\',
    'E:\scanner_data\twix_data\250221\',
    'E:\scanner_data\twix_data\250203\'};

MRFname{1}{1} = {'3D_MRF_fit_06-Mar-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD'};
MRFname{1}{2} = {'3D_MRF_DIP_07-Mar-2025_meas_MID00014_brain_HV5_FISP_FA15_nex600_w3000_3Dfit_SVD'};

MRFname{2}{1} = {'3D_MRF_fit_05-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD'};
MRFname{2}{2} = {'3D_MRF_DIP_07-Mar-2025_meas_MID00024_brain_HV4_FISP_FA15_nex600_3Dfit_SVD'};

MRFname{3}{1} = {'3D_MRF_fit_25-Feb-2025_meas_MID00547_brain_HV2_FISP_FA15_nex600_3Dfit_SVD'};
MRFname{3}{2} = {'3D_MRF_DIP_03-Mar-2025_meas_MID00547_brain_HV2_FISP_FA15_nex600_DIP_3Dvolume_HV2'};

MRFname{4}{1} = {'3D_MRF_fit_27-Feb-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_3Dfit_SVD'};
MRFname{4}{2} = {'3D_MRF_DIP_03-Mar-2025_meas_MID00561_brain_HV3_FISP_FA15_nex600_DIP_3Dvolume_HV3'};

MRFname{5}{1} = {'3D_MRF_fit_04-Feb-2025_meas_MID00188_brain_FISP_FA15_nex600_image_combined_3Dfit_SVD_tilted'};
MRFname{5}{2} = {'3D_MRF_DIP_12-Feb-2025_MID00188_DIP_3Dvolume'};

saveDir ='E:\POSTDOC_UoM\10_manuscripts\DIP_brainMRF_055T\data_HV_cluster\';
saveFolder = 'plots_cohort\';
mkdir([saveDir saveFolder])
save_filenames = 'plots_';

%%
for v = 1:length(MRFdir)
    for i = 1:2
        data{i} = load(strjoin([MRFdir{v} MRFname{v}{i} '.mat'],''));
        t1map(:,:,:,i,v) = squeeze(data{i}.t1map3d);
        t2map(:,:,:,i,v) = squeeze(data{i}.t2map3d);
        m0map(:,:,:,i,v) = squeeze(data{i}.m0map3d);
        b0map(:,:,:,i,v) = rot90(squeeze(data{i}.b0map3d));
    end
end


%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 250]; %uplimit:1200 for NIST T2 at 0.55T 

%% 1.1 test display
%v2 full
idx1{1} = 5:235; %vertical
idx2{1} = 65:255; %horizontal

idx1{2} = 5:235; %vertical
idx2{2} = 55:245; %horizontal

idx1{3} = 10:240; %vertical
idx2{3} = 55:245; %horizontal

idx1{4} = 5:235; %vertical
idx2{4} = 60:250; %horizontal

idx1{5} = 15:245; %vertical
idx2{5} = 40:230; %horizontal

toSliceViewer = squeeze(t1map(idx1{5},idx2{5},:,1,5));

% image = abs(vol3D(idx1,idx2,:));
% for n = 1:size(vol3D,3)
%     in = image(:,:,n);
%     image(:,:,n) = image(:,:,n)/max(in(:));
% end
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

%%
toplotLRT1 = cat(2,squeeze(t1map(idx1{1},idx2{1},:,1,1)),squeeze(t1map(idx1{2},idx2{2},:,1,2)),circshift(squeeze(t1map(idx1{3},idx2{3},:,1,3)),-11,3),circshift(squeeze(t1map(idx1{4},idx2{4},:,1,4)),-3,3),circshift(squeeze(t1map(idx1{5},idx2{5},:,1,5)),5,3));
toplotDIPT1 = cat(2,squeeze(t1map(idx1{1},idx2{1},:,2,1)),squeeze(t1map(idx1{2},idx2{2},:,2,2)),circshift(squeeze(t1map(idx1{3},idx2{3},:,2,3)),-11,3),circshift(squeeze(t1map(idx1{4},idx2{4},:,2,4)),-3,3),circshift(squeeze(t1map(idx1{5},idx2{5},:,2,5)),5,3));
toplotT1 = cat(1,toplotDIPT1,toplotLRT1);
figure,
sliceViewer(toplotT1, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar



%%
toplotLRT2 = cat(2,squeeze(t2map(idx1{1},idx2{1},:,1,1)),squeeze(t2map(idx1{2},idx2{2},:,1,2)),circshift(squeeze(t2map(idx1{3},idx2{3},:,1,3)),-11,3),circshift(squeeze(t2map(idx1{4},idx2{4},:,1,4)),-3,3),circshift(squeeze(t2map(idx1{5},idx2{5},:,1,5)),5,3));
toplotDIPT2 = cat(2,squeeze(t2map(idx1{1},idx2{1},:,2,1)),squeeze(t2map(idx1{2},idx2{2},:,2,2)),circshift(squeeze(t2map(idx1{3},idx2{3},:,2,3)),-11,3),circshift(squeeze(t2map(idx1{4},idx2{4},:,2,4)),-3,3),circshift(squeeze(t2map(idx1{5},idx2{5},:,2,5)),5,3));
toplotT2 = cat(1,toplotDIPT2,toplotLRT2);
figure,
sliceViewer(toplotT2, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar


%%
for v = 1:5
    c1m0 = abs(squeeze(m0map(:,:,:,1,v)));
    m0LR(:,:,:,v) = c1m0/max(c1m0(:))*1.8;
    c2m0 = abs(squeeze(m0map(:,:,:,2,v)));
    m0DIP(:,:,:,v) = c2m0/max(c2m0(:));
end
%%
toplotLRM0 = cat(2,m0LR(idx1{1},idx2{1},:,1),m0LR(idx1{2},idx2{2},:,2),circshift(m0LR(idx1{3},idx2{3},:,3),-11,3),circshift(m0LR(idx1{4},idx2{4},:,4),-3,3),circshift(m0LR(idx1{5},idx2{5},:,5),5,3));
toplotDIPM0 = cat(2,m0DIP(idx1{1},idx2{1},:,1),m0DIP(idx1{2},idx2{2},:,2),circshift(m0DIP(idx1{3},idx2{3},:,3),-11,3),circshift(m0DIP(idx1{4},idx2{4},:,4),-3,3),circshift(m0DIP(idx1{5},idx2{5},:,5),5,3));
toplotM0 = cat(1,toplotDIPM0,toplotLRM0);
figure,
sliceViewer(toplotM0, "DisplayRange", [0, 0.6], 'ScaleFactors', [1,1,1]), colorbar

%%
figure, histogram(toplotDIPM0)
th = 0.05;
for s = 1:size(c1m0,3)
    M0s = toplotDIPM0(:,:,s);
    cluster_masks(:,:,s) = M0s>th;
end

figure,
sliceViewer(cluster_masks), colorbar
%%
filename =[saveDir saveFolder save_filenames];
indexes = [11:37];
for n = 1:numel(indexes)
    slice = indexes(n);
    
    % Image = cat(4,squeeze(t1map(idx1,idx2,slice,:)),squeeze(t2map(idx1,idx2,slice,:)),squeeze(mm0(idx1,idx2,slice,:)));
    % ImageZ = cat(4,squeeze(t1map(idx1z,idx2z,slice,:)),squeeze(t2map(idx1z,idx2z,slice,:)),squeeze(mm0(idx1z,idx2z,slice,:)));
    
    fig=figure; 
    fig.Position = [0 0 800 1000];
    fig.Color = 'white';
    t = tiledlayout(3, 1, "TileSpacing", "none");
    
    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
    
    ax = nexttile;
    
    % ImageAux = Image(:,:,1,1);
    imagesc(toplotT1(:,:,slice).*cat(1,cluster_masks(:,:,slice),cluster_masks(:,:,slice)), climst1), colormap(ax, T1colormap)    
    axis off
    text(-0.03, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);    
    cb0 = colorbar(ax);
    cb0.Position(4) = cb0.Position(4)*0.9;
    cb0.FontSize = 10;

    ax = nexttile;
    imagesc(toplotT2(:,:,slice).*cat(1,cluster_masks(:,:,slice),cluster_masks(:,:,slice)), climst2), colormap(ax, T2colormap) 
    axis off
    text(-0.03, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    cb1.FontSize = 10;

    ax = nexttile;
    imagesc(toplotM0(:,:,slice).*cat(1,cluster_masks(:,:,slice),cluster_masks(:,:,slice)), [0, 0.6]), colormap(ax, 'gray')   
    axis off
    text(-0.03, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    cb2.FontSize = 10;

    exportgraphics(gcf,[saveDir saveFolder save_filenames '_s_' num2str(slice) '.tif'], 'Resolution',300)
    
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
