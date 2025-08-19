function [korig] = lowrank_2Dksp_inv(kcomp,Phi,idproj)
% Transforms low-rank MRF kspace back to original kspace
%
% Inputs
% 1. kcomp - (double, read x spiral_proj x coils x K) the low-rank MRF kspace,
% where K is the rank of the compressed MRF dictionary.
% 2. Phi -- (double, TRs x K) multiplication with this projection operator transforms data
% from the time domain to the compressed subspace. Multiplication by Phi'
% transforms data from the subspace back to the time domain.
% 4. idproj -- (double, TRs x 1) lists the spiral interleaf index for each TR
%
% Outputs
% 1. korig -- (double, read x coils x TRs) MRF k-space data
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[nr,np,nc,K]= size(kcomp);
Nex = size(Phi,1);
% imtx = Phi*Phi';
korig = complex(zeros(nr,nc,Nex,class(kcomp)));
for coil = 1:nc
    ksvd = reshape(kcomp(:,:,coil,:),[nr*np K])*Phi';
    for t=1:Nex
        p = idproj(t);
        korig(:,coil,t) = ksvd(nr*(p-1)+1:nr*p,t);
    end
end

