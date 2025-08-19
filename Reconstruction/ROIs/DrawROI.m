% Script to draw ROIs and get average and mean
% values on NIST phantom
% Jesus Fajardo jesuserf@umich.edu
filename = strcat('PM_SVDcompressed_Nex1000','.mat'); %mat filename
load(filename,'t1map3d','t2map3d','m0map3d');


t1map3d = ((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = ((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = ((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

slice = n;
hand_drawn_mask_phantom(t2map3d(:,:,slice), 14);
