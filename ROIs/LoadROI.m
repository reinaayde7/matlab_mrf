%{ 
Script to load Dicom files and hand-draw ROI 
and make statistics on that ROI for T1 and T2
jesuserf@med.umich.edu
%}

%% Load Dicom T1 and T2 values
filename1 = 'DICOM_T1map_B1Corr.dcm';
[info_sample1] = dicominfo(filename1);
[sample1] = dicomread(filename1);
T1vals = (squeeze(sample1(:,:,1,:)));
T1vals = double(T1vals);

filename2 = 'DICOM_T2map_B1Corr.dcm';
[info_sample2] = dicominfo(filename2);
[sample2] = dicomread(filename2);
T2vals = (squeeze(sample2(:,:,1,:)));
T2vals = double(T2vals);

load('PZ_ROI.mat', 'ROI');

%% Select slice
Slice = 16

data1 = T1vals(:,:,Slice).*ROI;
data2 = T2vals(:,:,Slice).*ROI;
data1 = nonzeros(data1);
data2 = nonzeros(data2);

format short
disp('mean T2: ')
sprintf('%16.f',mean(data2))
disp('Std T2: ')
sprintf('%16.f',std2(data2))

disp('mean T1: ')
sprintf('%16.f',mean(data1))
disp('Std T1: ')
sprintf('%16.f',std2(data1))

subplot()
climst1 = [0, 3000];
climst2 = [0, 800];

% Find boundary of logic matrix
boundary = bwboundaries(ROI);

imagesc(T2vals(:,:,Slice));%,climst1);
hold on
for k = 1:length(boundary)
    plot(boundary{k}(:,2), boundary{k}(:,1), 'r', 'LineWidth', 2)
end
hold off
colormap("gray");%(T1colormap);
title('ROI')
daspect([1 1 1]);
colorbar


%{
%imshow(histogram(data));
2800
%hand_drawn_mask_phantom(T2vals, 1);


%addpath('./Colormap');
figure('DefaultAxesFontSize',18)
subplot()
climst1 = [0, 3000];
climst2 = [0, 800];
imagesc(abs(quant_map));%,climst1);
colormap("gray");%(T1colormap);
title('M0')
colorbar

A = t2map3d(100:250,100:250,1);

subplot()
imagesc(A);
title('T2')
colorbar

A2 = 255*(A - min(A(:))) ./ (max(A(:)) - min(A(:))); %scale values between 0 and 255
A2 = cast(A2,'uint8');
imwrite(A2,'T22_deb.png','png'); 
%} 
