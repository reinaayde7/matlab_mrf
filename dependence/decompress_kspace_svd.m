function raw = decompress_kspace_svd(raw_svd, V, projID, fullview)
% Rudy 2501
% inputs
%   raw_svd: compressed data [ndata, fullview, nSVD, ncoils]
%   V: compression SVD matrix [Nframes, nSVD] (e.g., 600 7)
%   projID: idx of spiral proj (e.g. 0:47, 0:47, ...)
%   fullview: full number of projections per frame

[ndata, fullview, nk, ncoils] = size(raw_svd); %[4000 48 7 4] 
NFrames = size(V, 1); % [600]: number of acquired timepoints

raw_svd = permute(raw_svd, [3 1 2 4]); % [7 4000 48 4]
raw_svd = reshape(raw_svd, nk, []); % [7 768768]

% Invert the SVD compression
raw = ((raw_svd.')*V(:,1:nk)').'; %[768768 7]*[7 600] = [768768 600]' = [600 768768]
% raw = V(:,1:nk)*raw_svd; %[600 7]*[7 768768] = [600 768768]

% Reshape back to original dimensions
raw = reshape(raw, [NFrames, ndata, fullview, ncoils]); % [600 4000 48 4]

% Allocate space for the uncompressed raw data
raw_uncompressed = zeros(NFrames, ndata, ncoils); %[600,4000,4]

for t = 1:NFrames
    if isempty(projID)
        p = mod(t-1, fullview)+1;
    else
        p = projID(t)+1;
    end
    raw_uncompressed(t, :, :) = raw(t,:,p,:);
end

% Permute to match the original raw format
raw = permute(raw_uncompressed, [2, 1, 3]); % [4000,600,4]
end
