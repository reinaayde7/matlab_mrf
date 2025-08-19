tic
load("Compressed_Imgs_and_dict_Noisy3",'FirstComp','dict','r');
% Set the parameters for Gaussian noise

% Calculate the mean magnitude of the FirstComp array
mean_magnitude = mean(abs(FirstComp(:)));

% Define the percentage for the noise standard deviation
noise_percentage = .5; % For example, 10% of the mean magnitude
% Calculate the noise standard deviation as a percentage of the mean magnitude
stddev_noise = noise_percentage * mean_magnitude;
% Generate Gaussian noise for the real and imaginary parts
real_noise = stddev_noise * randn(size(FirstComp));
imaginary_noise = stddev_noise * randn(size(FirstComp));
% Add Gaussian noise to the complex array
noisy_FirstComp = FirstComp + (real_noise + 1i * imaginary_noise);
% Display the magnitude of the original and noisy complex arrays
figure;
subplot(1, 2, 1);
imshow((abs(FirstComp(:,:,3,16)) + 1), []); % Using log to enhance visibility
title('Original Complex Array Magnitude');
subplot(1, 2, 2);
imshow((abs(noisy_FirstComp(:,:,3,16)) + 1), []); % Using log to enhance visibility
title('Noisy Complex Array Magnitude');

[denoised,Sigma2,P,SNR_gain] = denoise_recursive_tensor(FirstComp,[5 5]);                                                                           %% RRRRRRR First Comp = [Nx, Ny, EIG, Nz] 
toc


t1map3d = zeros(size(FirstComp,1),size(FirstComp,1),size(FirstComp,4));
t2map3d = zeros(size(FirstComp,1),size(FirstComp,1),size(FirstComp,4));
m0map3d = zeros(size(FirstComp,1),size(FirstComp,1),size(FirstComp,4));
mask = true([size(FirstComp,1) size(FirstComp,1)]);
NumFrames = size(FirstComp,3);

Slice = 16;
[t1map_n,t2map_n,~,m0map_n] = patternmatch(noisy_FirstComp(:,:,:,Slice),mask,r,0,dict(1:NumFrames,:),4);
[t1map,t2map,~,m0map] = patternmatch(denoised(:,:,:,Slice),mask,r,0,dict(1:NumFrames,:),4);
save('DenoisedMaps_3.mat','t1map','t2map','m0map','-v7.3');
save('NoisyMaps.mat_3.mat','t1map_n','t2map_n','m0map_n','-v7.3');

quant_map = squeeze(rot90(abs(t1map_n(140:310,100:280))));
figure('DefaultAxesFontSize',16)
subplot()
climst1 = [0, 3000];
climst2 = [0, 1200];
climst3 = [-100, 100];
climst = [0, 1e-1];
imagesc(quant_map,climst1);
colormap(T1colormap);
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])
title('')
set(gca, 'position', [0 0 1 1], 'units', 'normalized')


quant_map = squeeze(rot90(abs(t1map(140:310,100:280))));
figure('DefaultAxesFontSize',16)
subplot()
climst1 = [0, 3000];
climst2 = [0, 1200];
climst3 = [-100, 100];
climst = [0, 1e-1];
imagesc(quant_map,climst1);
colormap(T1colormap);
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])
title('')
set(gca, 'position', [0 0 1 1], 'units', 'normalized')

quant_map = squeeze(rot90(abs(t2map_n(140:310,100:280))));
figure('DefaultAxesFontSize',16)
subplot()
climst1 = [0, 3000];
climst2 = [0, 1200];
climst3 = [-100, 100];
climst = [0, 1e-1];
imagesc(quant_map,climst2);
colormap(T2colormap);
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])
title('')
set(gca, 'position', [0 0 1 1], 'units', 'normalized')


quant_map = squeeze(rot90(abs(t2map(140:310,100:280))));
figure('DefaultAxesFontSize',16)
subplot()
climst1 = [0, 3000];
climst2 = [0, 1200];
climst3 = [-100, 100];
climst = [0, 1e-1];
imagesc(quant_map,climst2);
colormap(T2colormap);
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])
title('')
set(gca, 'position', [0 0 1 1], 'units', 'normalized')

%{
%quant_map = squeeze(rot90((t1map(130:280,160:310)-t1map_n(130:280,160:310))));%(160:310,115:265,Slice)(2*160:2*310,2*115:2*265,Slice)
quant_map = squeeze(rot90(abs(t2map(:,:))));%t2map(:,:)

figure('DefaultAxesFontSize',16)
subplot()
climst1 = [0, 3000];
climst2 = [0, 1200];
climst3 = [-100, 100];
climst = [0, 1e-1];
imagesc(quant_map,climst2);
colormap(T2colormap);
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])
title('')
set(gca, 'position', [0 0 1 1], 'units', 'normalized')


print(gcf, '-dpng', '-r1000', 'output.png'); % Optional: Save the figure if needed


save(tempfilename,'t1map3d','t2map3d','m0map3d','-v7.3');
subplot()
imagesc((abs(denoised(:,:,1,15))));
colormap("gray");
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])

subplot()
imagesc((abs(FirstComp(:,:,1,15))));
colormap("gray");
set(gca, 'XTick', [], 'YTick', []);
axis off;
daspect([1 1 1])

TR = 1000;
TI= 2000;
t1WeightingFunction = m0map.*( 1 - 2 * exp(-TI ./ t1map) + exp(-TR ./ t1map)); %this would be for an inversion-recovery
t1WeightedImage = t1map .* t1WeightingFunction;
t1WeightedImage = (t1WeightedImage - min(t1WeightedImage(:))) / (max(t1WeightedImage(:)) - min(t1WeightedImage(:)));

TE = 50;
t2WeightingFunction = m0map .* exp(-TE ./ t2map);
t2WeightedImage = t2map .* t2WeightingFunction;
t2WeightedImage = (t2WeightedImage - min(t2WeightedImage(:))) / (max(t2WeightedImage(:)) - min(t2WeightedImage(:)));
%}
