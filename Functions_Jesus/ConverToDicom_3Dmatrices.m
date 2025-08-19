filename = strcat('PM_SVDcompressed_SPcorrected_FSRdeblurr_Nex3000','.mat'); %mat filename
load(filename,'t1map3d','t2map3d','m0map3d');
%load('B1_map.mat','LPF_B1');
load('off_res.mat','LPF');

t1map3d = ((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = ((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = ((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));
%LPF_B1 = ((reshape(LPF_B1,[size(LPF_B1,1) size(LPF_B1,1) 1 size(LPF_B1,3)])));
LPF = ((reshape(LPF,[size(LPF,1) size(LPF,1) 1 size(LPF,3)])));

%m0map3d = flipud(m0map3d);

%hand_drawn_mask_phantom(t2map3d(:,:,1), 14);

%%%%%%%%%%%%%%DATA 2D MR7 from Customer Seq%%%%%%%%%%%%%%%%%%%%%%%% CHECKPOINT FEB 27 2023
%T1 = [2312.00, 1964.60, 1718.40, 1454.80, 1208.60, 960.40, 664.80, 517.40, 403.26, 270.80, 210.93, 164.86, 128.00, 92.71]
%T1SD = [76.02, 65.66, 106.53, 136.80, 110.18, 76.96, 42.58, 36.13, 15.21, 11.58, 7.81, 6.58, 6.08, 4.49]
%T2 = [579.20, 391.40, 253.00, 163.80, 117.60, 75.00, 51.44, 39.24, 25.52, 13.84, 7.64, 7.24, 4.36, 6.36]
%T2SD = [61.54, 27.18, 16.66, 11.14, 10.65, 7.99, 5.19, 3.70, 2.48, 4.73, 4.23, 3.30, 2.88, 4.53]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DATA 3D MR7 Data from Clinical prot.%%%%%%%%%%%%%%%%%%%%%%%% CHECKPOINT FEB 27 2023
%T1 = []
%T1SD = []
%T2 = []
%T2SD = []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data with NIST rotated 180º
%%%%%%%%%%%%%%DATA 2D MR7 from Customer Seq%%%%%%%%%%%%%%%%%%%%%%%% CHECKPOINT FEB 27 2023
%T1 = []
%T1SD = []
%T2 = []
%T2SD = []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DATA 3D MR7 Data from Clinical prot.%%%%%%%%%%%%%%%%%%%%%%%% CHECKPOINT FEB 27 2023
%T1 = []
%T1SD = []
%T2 = []
%T2SD = []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hand_drawn_mask_phantom(t2map3d(:,:,1), 14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file_name1 = strcat('DICOM_T1map','.dcm');
%file_name2 = strcat('DICOM_T2map','.dcm');
%file_name3 = strcat('DICOM_M0map','.dcm');

%dicomwrite(int16(t1map3d), file_name1);%
%dicomwrite(int16(t2map3d), file_name2); %(int16(t2map3d), file_name2, info, 'CreateMode', 'Copy', 'MultiframeSingleFile', 'true')
%dicomwrite(int16(m0map3d), file_name3);


%Crop image handrawn
%quant_map = squeeze(m0map3d(90:290,100:300,22));
%{
size_image_x = size(quant_map,1);
size_image_y = size(quant_map,2);
figure, imagesc(1:size_image_x, 1:size_image_y, quant_map);
h = drawfreehand;
ROI = createMask(h);
MaskedImage = ROI.*quant_map;
data = nonzeros(MaskedImage);
%histogram(data,15)
%imshow(MaskedImage,[]);
%hand_drawn_mask_phantom(t2map3d, 1);
%} 
% Define the number of slices to display
%{
D = t1map3d;
num_slices = 5;
climst1 = [0, 3000];
climst2 = [0, 800];
% Define the slice indices to display
slice_indices = [12, 13, 14, 15, 16, 17, 18];

% Create a figure to display the slices
figure;
% Set the spacing between subplots
% Set the horizontal spacing between subplots

% Loop over each slice index and plot the corresponding slice
for i = 1:length(slice_indices)
    % Get the slice index
    slice_index = slice_indices(i);
    
    % Extract the slice data
    slice_data = squeeze(D(150:300,100:200,slice_index));
    
    % Plot the slice data
    subplot(1, num_slices, i);
    imshow(slice_data, []);
    imagesc(rot90(slice_data,1),climst1);
    colormap(T1colormap);
    title('T1')
    daspect([1 1 1]) 
    %colorbar
    
    title(sprintf('Slice %d', slice_index));

end


%}
%addpath('./Colormap');
figure('DefaultAxesFontSize',2)
subplot()
climst1 = [0, 3000];
climst2 = [0, 800];
imagesc(rot90(squeeze(t2map3d(:,:,1))));%,1));%,climst2);
colormap(T2colormap);
title('T2')
colorbar

%{
A = t2map3d(100:250,100:250,1);

subplot()
imagesc(A);
title('T2')
colorbar

A2 = 255*(A - min(A(:))) ./ (max(A(:)) - min(A(:))); %scale values between 0 and 255
A2 = cast(A2,'uint8');
imwrite(A2,'T22_deb.png','png'); 
%}
