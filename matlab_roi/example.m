clear all, close all, clc
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))

%% step 1: get rois
% 1.1: start by plotting one of the coregistered maps
load('MRFMaps_DirectMatchSVD.mat');
figure;
imagesc(t1SVD);
axis image;

% 1.2: zoom in to the relevant region, make sure you unselect the 
% zoom option when you're done!

%% 1.3: draw rois
% execute the following command. Then go to your figure again. Type 'c' 
% to draw a circular ROI, 'f' to draw a freehand ROI. Each new ROI will be 
% converted to a mask and added to a cell array. Type 'q' when you are 
% done. If you later want to compare to NIST values, make sure you start at
% the highest value (brightest contrcast sphere). 
rois = get_rois;

%% step 2: plot against NISTq
NIST = load('NIST_T2layer_1,5T.txt');
plot_vs_NIST(NIST, t1SVD, t2SVD, rois);

%% Alternative: investigate ROIS individually
% You can also check values in an ROI individually. Works very similar to 
% drawing ROIs. Just give any coregistered arrays you want to investigate
% as arguments and you'll get means and stds immediately in the terminal. 
% Draw ROIs as you did in step 1.3 after zooming in on the relevant region.
% To exit, just close the figure.
figure;
imagesc(t1SVD);
axis image;
