%{ 
Script to load Dicom files and hand-draw ROI from one of them, save
the mask and use it in a second array
Jesus Fajardo jesuserf@umich.edu
%}
%% Load Dicom T1 and T2 values
filename1 = 'DICOM_T1map.dcm';
[info_sample1] = dicominfo(filename1);
[sample1] = dicomread(filename1);
T1vals = (squeeze(sample1(:,:,1,:)));

filename2 = 'DICOM_T2map.dcm';
[info_sample2] = dicominfo(filename2);
[sample2] = dicomread(filename2);
T2vals = (squeeze(sample2(:,:,1,:)));

%% Select slice
Slice = 16

%Crop image handrawn
quant_map = squeeze(T2vals(:,:,Slice));
%hand_drawn_mask_phantom(quant_map, 1);

size_image_x = size(quant_map,1);
size_image_y = size(quant_map,2);
figure, imagesc(1:size_image_x, 1:size_image_y, quant_map);
h = drawfreehand;
ROI = createMask(h);
MaskedImage2 = int16(ROI).*quant_map;
data2 = nonzeros(MaskedImage2);
data1 = nonzeros(int16(ROI).*T1vals(:,:,Slice));
imshow(ROI);%(int16(ROI).*T1vals(:,:,Slice),[])

save('PZ_ROI.mat', 'ROI');

format short
disp('mean T2: ')
sprintf('%16.f',mean(data2))
disp('Std T2: ')
sprintf('%16.f',std2(data2))

disp('mean T1: ')
sprintf('%16.f',mean(data1))
disp('Std T1: ')
sprintf('%16.f',std2(data1))
