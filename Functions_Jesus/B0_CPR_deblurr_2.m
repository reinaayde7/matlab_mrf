function final_image = B0_CPR_deblur_2(images, demod_freqs, fieldmap);
% PART TWO (II) of CPR spiral image deblurring.
% Takes images reconstructed at L+1 freqs 
% and the coefficients from B0_CPR_deblur_1 and
% reconstructs the final deblurred image
% with the pixel of the demodulated image
% closest to the B0 image value
% Jesus Fajardo (jesuserf@med.umich.edu)
disp('CPR-Deblurring...')

demod_freqs = demod_freqs / (2*pi);
field_map = fieldmap;
final_image = zeros(size(field_map,1), size(field_map,2),size(images,3));

for i = 1:size(field_map,1)
    for j = 1:size(field_map,2)
        [ d, idx ] = min( abs(demod_freqs - field_map(i,j)));

        %disp('DeltaB0:')
        %disp(field_map(i,j))
        %disp('Freqs:')
        %disp(demod_freqs)
        %disp('Actual value:')
        %disp(demod_freqs(idx))
        %disp('index found:')
        %disp(idx)
        final_image(i,j,:) = images(i,j,:,idx);
    end
end


end