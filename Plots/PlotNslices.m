% Script to plot N slices from 3D matrices
% Jesus Fajardo jesuserf@umich.edu

filename = strcat('PM_SVDcompressed_Nex1000','.mat'); %mat filename
load(filename,'t1map3d','t2map3d','m0map3d');
t1map3d = rot90((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = rot90((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = rot90((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

addpath('./Colormap');

D = (t2pf(:,:,:));
montage(D, "Indices",40:45,"DisplayRange",[0 300], 'hot');


