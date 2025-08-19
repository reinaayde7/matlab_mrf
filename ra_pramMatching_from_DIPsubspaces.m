% Reina Ayde
% pattern matching giving DIP subspaces of one slice
%
% Based on previous scripts by Rudy Rizzo
% Based on previous scripts by Jesus Fajardo (jesuserf@umich.edu)
% Based on previous scripts by Yun Jiang (yunjiang@med.umich.edu)
clear;clc;
% tic
%% setup path and irt

addpath('./rr_dictionary/DictSimulation_NoSP');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');
% import npy-matlab package
addpath('C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\npy-matlab\npy-matlab');


% Path to IRT toolbox
irtPath = 'C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\irt\irt';
addpath(irtPath); run setup;

codePath = 'C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\imagine\'; % ra: I commented momentarily
addpath(genpath(codePath)) % ra: I commented momentarily

%% if import the subspaces: numpy file of IRN_output
%path = 'C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\matlab_scripts\3D_MRF_FISP_Prostate-main_d240531-main\results\recon_25-Apr-2025_trueFISP\DIP\FISP\IRN_output_s5.npy'
path = 'C:\Users\ayde\PycharmProjects\dip_tensorflow\results_from_cluster\20250521_IRNoutputs_NIST\'
%image_combined_double = load(fullfile(path)).subspaces; %(304 x 304 x RANK)
%image_combined = squeeze(image_combined_double(3, :, :, :));
IRN_output = 'IRN_output_epoch2.npy'
image_combined_double = readNPY(fullfile(path, IRN_output));
image_combined = complex(image_combined_double(:, :, 1:ceil(size(image_combined_double, 3)/2)), image_combined_double(:, :, ceil(size(image_combined_double, 3)/2)+1:end));
%% Load corresponding dictionary
dict = load(fullfile(path, 'dictCompressed.mat'));
dictCompressed = dict.dictCompressed;
Phi = dict.Phi;
Vc = Phi;
r = dict.r;
NumFrames = size(dictCompressed, 1);
%% Pattern matching low rank data
fprintf('LR recon\n');
FOV = 304 ;
mask_fit = true([FOV FOV]);
tic
[t1map,t2map,b0map,m0map] = patternmatch(image_combined, mask_fit,r,0,dictCompressed(1:NumFrames,:),16); 
toc

%% single slice plots
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500]; %3000 NIST, 1500 brain, 2500 heart
climst2 = [0, 300];  %1500 NIST, 300 brain, 150 heart
% idx1 = [100:200];
% idx2 = [100:200];
s_vert = 40:260;
s_hor = 60:300;
idx1 = s_vert;
idx2 = s_hor;
figure, 
ax1 = subplot(311); imagesc(rot90(t1map(idx1,idx2)), climst1), colormap(ax1, T1colormap), colorbar
ax2 = subplot(312); imagesc(rot90(t2map(idx1,idx2)), climst2), colormap(ax2, T2colormap), colorbar

m00 = abs(m0map)/max(abs(m0map(:)));
%m00LR = abs(m0mapLLR)/max(abs(m0mapLLR(:)));
ax3 = subplot(313); imagesc(rot90(m00(idx1,idx2)),[0 0.6]), colormap(ax3, gray), colorbar
