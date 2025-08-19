function dixon_formatted = B0_CPR_deblur_1(raw_arr);
% Formats vibe_dixon from Siemens
% to the input required by HierarchicalIDEAL
% Nx Ny Nslices Ncoils NTE
% Jesus Fajardo (jesuserf@med.umich.edu)
%set(gca,'DefaultTextFontSize',28);
%clear;clc;close all;
cd('..');
addpath('./dependence');
addpath('./OpenSiemensRawData');
addpath('./Functions_Jesus/');
addpath('./data/');
%addpath('./data/896/');
%cd('..');


%% Open Siemens 2D Raw data
filename = 'data/Dixon_6echo_100Hz'
load('Liver_ax.mat','r')    

raw = mapVBVD(filename);
raw = raw{2}; %not using noise data

raw = squeeze(raw.image( :, :, :, :, :, :, :, :, :, :, :)); % Select slice in 5th pos
raw = permute(raw,[1,3,4,2,5]);

data = struct('images',raw, ...
              'TE', [0.00093 0.00235 0.0037 00519 00661 00803], ... 
              'FieldStrength', 1.5000, ...
              'PrecessionIsClockwise', 1);

save('Phantom_Aug_23_6Echo_100hz_2.mat','data')





end