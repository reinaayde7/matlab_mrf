function data = LoaDixon(filename)%, maskbase)

% function LoaDixon = LoaDixon( filename )
% Script to load Fat and water maps using
% dixon sequence.
% filename is the .dat raw Siemens data
% Jesus Fajardo (jesuserf@med.umich.edu)

set(gca,'DefaultTextFontSize',28);
%clear;clc;close all;
cd('..');
addpath('./dependence');
addpath('./OpenSiemensRawData');
addpath('./Functions_Jesus/');
addpath('./data/');
%addpath('./data/896/');
%cd('..');

%% Open Siemens 2D Raw data
filename = 'data/meas_MID03186_FID126966_T1_vibe_Dixon_FOV400'

raw = mapVBVD(filename);
raw = raw{2}; %not using noise data

%disp(size(raw))

%raw.image

raw = squeeze(raw.image( :, :, :, :, :, :, :, :, :)); % Select slice in 5th pos
data = permute(raw,[3,1,4,2,5]);
[Nx,Ny,Nslices,Nc,NTE] = size(data);


disp('raw data size:')
disp(size(raw))


%% Classic FFT
disp('Do fft through-plane')
data = transformKspaceToImage(data,3);
fprintf('FFT...\n')
Nslice = floor(Nslices / 2); %Here I'm selecting just a single slice, otherwise it runs outta memory
%A single slice imageNTE [Nx, Ny, Nc]:
image1 = zeros(size(data,1),size(data,2),size(data,4));
image2 = zeros(size(data,1),size(data,2),size(data,4));
image3 = zeros(size(data,1),size(data,2),size(data,4));
image4 = zeros(size(data,1),size(data,2),size(data,4));
image5 = zeros(size(data,1),size(data,2),size(data,4));
image6 = zeros(size(data,1),size(data,2),size(data,4));


for i = 1:Nc
    image1(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,1))));
    image2(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,2))));
    image3(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,3))));
    image4(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,4))));
    image5(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,5))));
    image6(:,:,i) = fftshift(ifft2(squeeze(data(:,:,Nslice,i,6))));
end

%% Coil Sensitivity Map estimation
fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

avg_image1 = sum(image1,3);
avg_image2 = sum(image2,3);
avg_image3 = sum(image3,3);
avg_image4 = sum(image4,3);
avg_image5 = sum(image5,3);
avg_image6 = sum(image6,3);

csm1 = estimate_csm_walsh(avg_image1);
csm2 = estimate_csm_walsh(avg_image2);
csm3 = estimate_csm_walsh(avg_image3);
csm4 = estimate_csm_walsh(avg_image4);
csm5 = estimate_csm_walsh(avg_image5);
csm6 = estimate_csm_walsh(avg_image6);

%subplot();
%imagesc(abs(squeeze(csm1(:,:,8))));

%% Perform coil combination
image_combined1 = single(zeros(size(image1,1),size(image1,3),Nc));
image_combined2 = single(zeros(size(image2,1),size(image2,3),Nc));
image_combined3 = single(zeros(size(image3,1),size(image3,3),Nc));
image_combined4 = single(zeros(size(image4,1),size(image4,3),Nc));
image_combined5 = single(zeros(size(image5,1),size(image5,3),Nc));
image_combined6 = single(zeros(size(image6,1),size(image6,3),Nc));

image_combined1 = sum(image1.*conj(csm1),3)./sum(csm1.*conj(csm1),3); 
image_combined2 = sum(image2.*conj(csm2),3)./sum(csm2.*conj(csm2),3); 
image_combined3 = sum(image3.*conj(csm2),3)./sum(csm3.*conj(csm3),3); 
image_combined4 = sum(image4.*conj(csm2),3)./sum(csm4.*conj(csm4),3); 
image_combined5 = sum(image5.*conj(csm2),3)./sum(csm5.*conj(csm5),3); 
image_combined6 = sum(image6.*conj(csm2),3)./sum(csm6.*conj(csm6),3); 

image_combined1 = image_combined1(:,size(image_combined1,2)/4:3*size(image_combined1,2)/4-1); % center image and (JESUS added -1)
image_combined2 = image_combined2(:,size(image_combined2,2)/4:3*size(image_combined2,2)/4-1); % remove oversampling
image_combined3 = image_combined3(:,size(image_combined3,2)/4:3*size(image_combined3,2)/4-1);
image_combined4 = image_combined4(:,size(image_combined4,2)/4:3*size(image_combined4,2)/4-1);
image_combined5 = image_combined5(:,size(image_combined5,2)/4:3*size(image_combined5,2)/4-1);
image_combined6 = image_combined6(:,size(image_combined6,2)/4:3*size(image_combined6,2)/4-1);

inputimg = zeros(size(image_combined1,1),size(image_combined1,2),1,1,NTE);

inputimg(:,:,1,1,1) = image_combined1;
inputimg(:,:,1,1,2) = image_combined2;
inputimg(:,:,1,1,3) = image_combined3;
inputimg(:,:,1,1,4) = image_combined4;
inputimg(:,:,1,1,5) = image_combined5;
inputimg(:,:,1,1,6) = image_combined6;


subplot()
imagesc(abs(squeeze(inputimg(:,:,1,1,1))));
title('Magnitude Echo1')
colorbar

data = struct('images',inputimg, ...
              'TE', [0.0091 0.00235 0.00374 0.0513 0.0652 0.0791], ... 
              'FieldStrength', 3.00, ...
              'PrecessionIsClockwise', 1);

save('Phantom_Oct_24_6Echo_300hz.mat','data')

end
