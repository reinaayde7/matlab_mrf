function [demod_data, coeff_table] = B0_MFI_deblur_1(vector, field_map, L);
% PART ONE (I) of spiral image deblurring.
% Takes raw k-space spiral data, field map, L and 
% trajectory and returns demodulated raw k-space
% data and MFI c_i coefficients.
% Jesus Fajardo (jesuserf@med.umich.edu)


% Here vector is a raw array with dimensions [LenghtReadout,Nshots,NumCoils,NumFrames]

[LenghtReadout,Nshots,NumCoils,NumFrames] = size(vector);
ts = 2.5e-6; % time step

disp('read time')
read_time = LenghtReadout * ts;

%% DEFINE FREQUENCY SEGMENTS
%======================================================================
delta = (max(field_map(:)) - min(field_map(:)))/L*2*pi; % freq spacing
deltaw = linspace(max(field_map(:)), min(field_map(:)),L+1)*2*pi; % array of frequencies
[h,w] = size(field_map);

%% DEMODULATE DATA AT FREQUENCIES PRIOR TO GRIDDING
%======================================================================
demod_data = zeros(LenghtReadout,Nshots,NumCoils,NumFrames,L+1);
time = linspace(0,read_time,LenghtReadout);

disp('demodulated omegas')
disp(deltaw)
for ii = 1:L+1
    demod_data(:,:,:,:,ii) = vector.* exp(-1i * deltaw(ii) * time(:));
    disp('mean demod data number')
    disp(ii)
    disp(mean(demod_data(:,:,:,:,ii), 'all'))
end

   
%% GENERATE INTERPOLATION COEFFICIENTS
%========================================
num_of_interp = 1000;%round(1000/L);
tk = 64;
timesamples = linspace(0,read_time,tk);
interpolation_omegas = linspace(min(field_map(:)),max(field_map(:)),num_of_interp)*2*pi;
% Approximate the linear combination coefficients using least squares
deltawi_tk = exp(1i * timesamples.' * deltaw );
y = exp(1i*timesamples.'*interpolation_omegas);
coeff_table = (deltawi_tk \ y).';
%disp(size(coeff_table));

end