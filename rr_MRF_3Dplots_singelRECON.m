% rr MRF 3D plot and gifs
% ONLY FOR ONE FILE!!
%Rudy, December 2024

clear all, close all, clc

doDebug = 0;
doSave = 0;
    savedir = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\';
    name = 'carotid_volume1_LR.gif';
%%

%volunteer 241018 3D LR recon 

MRFdir = ['E:\POSTDOC_UoM\08_Project_MRF\DIP\FISP\'];
MRFname = ['3D_MRF_DIP_10-Jan-2025_carotid_vol2_FISP_Nex1000'];

% MRFdir = ['E:\scanner_data\twix_data\241217\'];
% MRFname = ['3D_MRF_fit_18-Dec-2024_meas_MID00043_carotid_FISP_FA15_nex1000_np18'];


data = load([MRFdir MRFname '.mat']);
t1map = squeeze(data.t1map3d);
t2map = squeeze(data.t2map3d);
m0map = squeeze(data.m0map3d);
b0map = rot90(squeeze(data.b0map3d));

%%
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 200]; %uplimit:1200 for NIST T2 at 0.55T 

%% 1.1 test display

%volume 1 full
idx1 = 50:230; %vertical
idx2 = 60:260; %horizontal

% zoom in
idx1z = 120:160; %vertical
idx2z = 180:220; %horizontal

%volume 2 full
% idx1 = 30:220; %vertical
% idx2 = 60:260; %horizontal

% zoom in
% idx1z = 110:150; %vertical
% idx2z = 175:215; %horizontal

toSliceViewer = squeeze(t1map(idx1,idx2,:));

figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar
%% 1.2 video generator
% design from ISMRM25 where 3 different reconstructions where compared

filename =[savedir name];

for n = 1:8 % keep this loop running from 1 even though you plan to save other partition (otherwise error from gif generator)
    slice = n+5;
    mm0 = squeeze(abs(m0map(:,:,slice)));
    M = max(mm0(:));
    mm0 = mm0./M;
    
    % figure,
    Image = cat(3,squeeze(t1map(idx1,idx2,slice)),squeeze(t2map(idx1,idx2,slice)),mm0(idx1,idx2));
    ImageZ = cat(3,squeeze(t1map(idx1z,idx2z,slice)),squeeze(t2map(idx1z,idx2z,slice)),mm0(idx1z,idx2z));
    
    fig=figure; 
    fig.Position = [0 0 800 800];
    t = tiledlayout(3, 2, "TileSpacing", "tight");
    
    %--------------------------------------------------------------------
    % ROW 1

    %T1 maps
    colormapP = T1colormap;
    climst = climst1;
    
    %---------------------------------
    % COL 1
    ax = nexttile;  
    ImageAux = Image(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % axis equal
    pbaspect([1.5 1 1])
    axis off
    text(-0.1, 0.5, 'T1 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    % text(0.5, 1.1, 'LR-FISP', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 0);
    rectangle('Position', [115, 75, 45, 50], 'EdgeColor', 'y', 'LineWidth', 2);
    
    %---------------------------------
    % COL 2
    ax = nexttile;  
    ImageAux = ImageZ(:,:,1);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    % axis equal
    pbaspect([1 0.9 1])
    axis off

    cb0 = colorbar(ax);
    % cb0.Position(1) = cb0.Position(1)+0.05;
    cb0.Position(4) = cb0.Position(4)*0.9;
    
    %--------------------------------------------------------------------
    % ROW 2

    %T2 maps
    colormapP = T2colormap;
    climst = climst2;
    
    %---------------------------------
    % COL 1
    ax = nexttile;
    ImageAux = Image(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    pbaspect([1.5 1 1])
    axis off
    text(-0.1, 0.5, 'T2 [ms]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    
    %---------------------------------
    % COL 2
    ax = nexttile;
    ImageAux = ImageZ(:,:,2);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    pbaspect([1 0.9 1])
    axis off
    
    cb1 = colorbar(ax);
    cb1.Position(4) = cb1.Position(4)*0.9;
    
    %--------------------------------------------------------------------
    % ROW 3

    % M0 (or rho)
    colormapP = 'gray';
    climst = [0,1];

    %---------------------------------
    % COL 1
    ax = nexttile;
    ImageAux = Image(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    pbaspect([1.5 1 1])
    axis off
    text(-0.1, 0.5, 'M0 [a.u.]', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
    
    %---------------------------------
    % COL 2
    ax = nexttile;
    ImageAux = ImageZ(:,:,3);
    imagesc(ImageAux, climst), colormap(ax, colormapP)
    pbaspect([1 0.9 1])
    axis off
   
    cb2 = colorbar(ax);
    cb2.Position(4) = cb2.Position(4)*0.9;
    
    
    if doDebug
        break;
    end
    if doSave
        %gif saves
        timer=1;
        frame = getframe(fig);
        im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if n == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',timer);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',timer);
            end
    end
end
