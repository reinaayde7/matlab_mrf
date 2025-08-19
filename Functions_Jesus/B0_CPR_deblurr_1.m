function [demod_data, deltaw] = B0_CPR_deblur_1(vector, field_map, L);
% PART ONE (I) of spiral image deblurring.
% Takes raw k-space spiral data, field map, L and 
% trajectory and returns demodulated raw k-space
% data and MFI c_i coefficients.
% Jesus Fajardo (jesuserf@med.umich.edu)


% Here vector is a raw array with dimensions
% [LenghtReadout,Nshots,NumCoils,NumFrames] %Pay attention if you're doing
% SVD first and change this dimensions

[LenghtReadout,Nshots,NumCoils,NumFrames] = size(vector);
ts = 2.5e-6; % time step
read_time = LenghtReadout * ts;% * 1.e-3;

%% DEFINE FREQUENCY SEGMENTS
%======================================================================
delta = (max(field_map(:)) - min(field_map(:)))/L*2*pi; % freq spacing
deltaw = linspace(max(field_map(:)), min(field_map(:)),L+1)*2*pi; % array of frequencies
[h,w] = size(field_map);

%{
if L < floor(8*max(deltaw)*read_time)
    disp('The number of demodulating frequencies needs to be higher than')
    disp(floor(8*max(deltaw)*read_time))
    return
end
%}
%% DEMODULATE DATA AT FREQUENCIES PRIOR TO GRIDDING
%======================================================================
demod_data = zeros(LenghtReadout,Nshots,NumCoils,NumFrames,L+1);
time = linspace(0,read_time,LenghtReadout);

disp('Demodulating data at frequencies...')
for ii = 1:L+1
    demod_data(:,:,:,:,ii) = vector.* exp(-1i * deltaw(ii) * fliplr(time(:)));
end


end