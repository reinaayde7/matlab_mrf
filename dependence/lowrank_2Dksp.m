function [kcomp ksvd] = lowrank_2Dksp(korig,idproj,Phi,np)
% Compresses MRF k-space along the time dimension,
% then grids the compressed k-space to the image domain. Yields "low-rank"
% kspace
% 
% Input
% 1. korig -- (double, read x coils x TRs) the MRF k-space data
% 3. idproj -- (double, TRs x 1) spiral interleaf index for each TR
% 4. Phi -- (double, TRs x K) projection matrix that converts data from the
%               time domain to the compressed subspace. K is the rank of
%               the dictionary after compression.
% 6. np -- (double, 1x1) the number of spiral interleaves in the fully sampled data
% 
% Outputs
% kcomp -- (double, read x spiral_proj x coils x K) the low-rank kspace 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[nr,nc,Nex] = size(korig);
K = size(Phi,2);

% Compress k-space along the time dimension
ksvd = complex(zeros(nr*np,Nex,class(korig)));
kcomp = complex(zeros(nr,np,nc,K,class(korig)));
for coil = 1:nc
    for t=1:Nex
        p = idproj(t); % spiral interleaf
        ksvd(nr*(p-1)+1:nr*p,t) = korig(:,coil,t);
    end
    kcomp(:,:,coil,:) = reshape(ksvd*Phi,[nr np K]);
end

