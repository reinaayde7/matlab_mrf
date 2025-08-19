function raw_svd = compress_kspace_svd(raw,V,projID,fullview)
% inputs
%   raw: data [ndata, Nframes, ncoils] (e.g., 4000, 600, 4]
%   V: compression svd matrix [Nframes, nSVD] (e.g., 600 7)
%   projID: idx of spiral proj (e.g. 0:47, 0:47, ...)
% fullview = 48;

if ndims(raw)==3 %single partition or 2D slice
    [ndata,NFrames,ncoils] = size(raw); % [4000,600,4]
    nk = size(V,2); % [7]: number of SDV subspaces
       
    raw = permute(raw,[2 1 3]); %[600,4000,4]
    raw_svd = zeros(NFrames,ndata,fullview,ncoils); % [600 4000 48 4]
    for t=1:NFrames %sorting data by spiral projection (ie, zerofilling on kspace)
        
        if isempty(projID)
            p = mod(t-1,fullview)+1;
        else
            p = projID(t)+1;
        end
        raw_svd(t,:,p,:) = raw(t,:,:);
    end; clear raw;

    raw_svd = reshape(raw_svd,NFrames,[]); % [600 768768]
    %disp(size(raw_svd')) JESUS
    %disp(size(V(:,1:nk))) JESUS

    % debug.r1 = raw_svd; %[600 768768]

    % raw_svd = ((raw_svd.')*V(:,1:nk)).'; %[768768 600]*[600 7] = [768768 7]' = [7 768768]
    
    raw_svd = (V(:,1:nk).')*raw_svd; %[7 600]*[600 768768] = [7 768768]
    
    % debug.r2 = raw_svd;
    % debug.r3 = V(:,1:nk)'*(V(:,1:nk)*raw_svd); %[7 600] [600 7] [7 768768] = [600 768768]

    % raw = (raw*V(:,1:nk)).';
    raw_svd = reshape(raw_svd,[nk ndata fullview ncoils]); %[7 4000 48 4]
    raw_svd = permute(raw_svd,[2 3 1 4]);%[4000 48 7 4] 
    
elseif ndims(raw)==4
    [ndata,nz,NFrames,ncoils] = size(raw);
    nk = size(V,2);
    
    raw = permute(raw,[3 1 2 4]);
    raw_svd = zeros(NFrames,ndata,fullview,nz,ncoils);
    for t=1:NFrames
        if isempty(projID)
            p = mod(t-1,fullview)+1;
        else
            p = projID(t);
        end; 
    end;clear raw;
    raw_svd = reshape(raw_svd,NFrames,[]);
    raw_svd = ((raw_svd.')*V(:,1:nk)).';
    % raw = (raw*V(:,1:nk)).';
    raw_svd = reshape(raw_svd,[nk ndata ncoils,nz]);
    raw_svd = permute(raw_svd,[2 3 4 1]);
% elseif ndims(raw)==5
%     [ndata,nshots,nz,NFrames,ncoils] = size(raw);
%     nk = size(V,2);
%     
%     raw = permute(raw,[4 1 2 3,5]);
%     raw_svd = zeros(NFrames,ndata,nshots,fullview,ncoils,nz);
%     for t=1:NFrames
%         p = projID(t);
%         raw_svd(t,:,p,:,:) = raw(t,:,:,:);
%     end; clear raw;
%     raw_svd = reshape(raw_svd,NFrames,[]);
%     raw_svd = ((raw_svd.')*V(:,1:nk)).';
%     % raw = (raw*V(:,1:nk)).';
%     raw_svd = reshape(raw_svd,[nk ndata ncoils,nz]);
%     raw_svd = permute(raw_svd,[2 3 4 1]);
else
    error('image matrix has wrong size')
end
