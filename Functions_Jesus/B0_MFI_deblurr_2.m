function final_image = B0_MFI_deblur_2(images, coeff, fieldmap, L, num_of_interp);
% PART TWO (II) of MFI spiral image deblurring.
% Takes images reconstructed at L+1 freqs 
% and the coefficients from B0_MFI_deblur_1 and
% reconstructs the final deblurred image
% as linear combination
% Jesus Fajardo (jesuserf@med.umich.edu)


%Test parameters
%L = 5;
%load('off_res.mat','masked3'); % B0 map (off-res freq)
%load('c_i.mat','coeff_table'); % Coefficients
%load('Demod_Imgs.mat','image_uncombined'); % Demodulated Images
coeff_table = coeff;
field_map = fieldmap;
interpolation_omegas = linspace(min(field_map(:)),max(field_map(:)),num_of_interp)*2*pi;
%images = image_uncombined;

%% Create final image
%========================================
%Add the pixels together over the images using the interpolation
%coefficients
[fh, fw] = size(field_map);
%permute images for faster calcualation
%images = permute(images,[2 3 1]);
final_image = zeros(fh,fw);
[~,~,Ncoils,Nshots,~] = size(images);
final_image = zeros(fh,fw,Ncoils,Nshots);


%pause(99999999999999999)
disp('reconstructing final MFI corrected image...')
%for nco = 1:Ncoils
    for ns = 1: Nshots
        for x = 1:fw
                for y = 1:fh
                    px_omega = field_map(x,y)*2*pi;
                    [~, idx] = min(abs(interpolation_omegas - px_omega));
                    final_image(x,y,ns) = coeff_table(idx,:) * squeeze(images(x,y,:,ns));
                    %index_progress = Ncoils+Nshots+fw+fh;
                end
            
        end
    end
%end


%% Normalisation of the image
mi = min(final_image(:));
final_image = final_image - mi;
ma = max(final_image(:));
final_image = final_image/ma;

end