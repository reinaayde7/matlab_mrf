% Coil PCA compression
%
% V0.1 - Modified based on Jesse's original code
%      - kdata has a fixed number of dimension.
%      - arrange as [nADCpoints,nshots,nframes,ncoils] for 2D dataset
%      - 12.06.2015
% v0.2 - modified to accodemate prostate multi-shot data -- 01/25/2016
% Yun Jiang -- yun.jiang@case.edu

function [kdata,nVirtualCoils,V] = pcaCoilCompress(kdata,thresh,varargin)
% kdata: assumed to be single-shot 2D
% thresh: a number between 0-1 that defines what percentage of the original coil variations are retained
% 0.95 is a good choice, and will keep 95% of the coil sensitivity variations
%
% OUTPUTS:
% kdata: the compressed k-space matrix [nADCpoints,ncoils,nshots,nframes]
% nVirtualCoils: the number of coils after compression
% V: eigenvector matrix (save this if you need to compress other datasets the same way)
if ndims(kdata) == 4,
    [nADCpoints,nshots,nframes,ncoils] = size(kdata); %
    %[nADCpoints,nshots,nframes,ncoils] = size(kdata);
    D = reshape(kdata,[nADCpoints*nshots*nframes],ncoils);
    
    if isempty(varargin)
        
        [U,S,V] = svd(D,'econ');
        coilinfo = cumsum(diag(S))/sum(diag(S));
        
        %         figure(1);clf; plot(1:nc,coilinfo); xlabel('Coils'); ylabel('SVD'); title('PCA Coil Compression');
        nVirtualCoils = find(coilinfo > thresh,1,'first');
        fprintf('\t Using %d virtual channels (%d coils originally)\n',nVirtualCoils,ncoils);
        
    elseif length(varargin)==2
        V = varargin{1};
        nVirtualCoils = varargin{2};
    else
        nVirtualCoils = varargin{1};
        [U,S,V] = svd(D,'econ');
    end
    
    temp = zeros(nADCpoints,nshots,nframes,nVirtualCoils,'single');
    for i=1:nframes
        for j = 1:nshots
            temp(:,j,i,:) = reshape(reshape(kdata(:,j,i,:),nADCpoints,ncoils)*V(:,1:nVirtualCoils),[nADCpoints nVirtualCoils]);
        end
    end
    kdata = temp;clear temp;
    
elseif ndims(kdata) == 5,
    [nADCpoints,nshots,npartitions,nframes,ncoils] = size(kdata); %
    %[nADCpoints,nshots,nframes,ncoils] = size(kdata);
    D = reshape(kdata,[nADCpoints*nshots*npartitions*nframes],ncoils);
    
    if isempty(varargin)
        
        [U,S,V] = svd(D,'econ');
        coilinfo = cumsum(diag(S))/sum(diag(S));
        
        %         figure(1);clf; plot(1:nc,coilinfo); xlabel('Coils'); ylabel('SVD'); title('PCA Coil Compression');
        nVirtualCoils = find(coilinfo > thresh,1,'first');
        fprintf('\t Using %d virtual channels (%d coils originally)\n',nVirtualCoils,ncoils);
        
    elseif length(varargin)==2
        V = varargin{1};
        nVirtualCoils = varargin{2};
    else
        nVirtualCoils = varargin{1};
        [U,S,V] = svd(D,'econ');
    end
    
    temp = zeros(nADCpoints,nshots,npartitions,nframes,nVirtualCoils,'single');
    for i=1:nframes
        for j = 1:nshots
            for ipar = 1:npartitions
                temp(:,j,ipar,i,:) = reshape(reshape(kdata(:,j,ipar,i,:),nADCpoints,ncoils)*V(:,1:nVirtualCoils),[nADCpoints nVirtualCoils]);
            end
        end
    end
    kdata = temp;clear temp;
end;


% elseif ndims(kdata)==5
%
%     %% pca coil compression
%     fprintf('PCA coil compression...\n')
%     % kdata originally has size: read x proj x coils x part x meas
%     kdata = permute(kdata,[1 2 4 5 3]);
%     [nr,npacc,nz,nttot,nc] = size(kdata);
%
%     %     kdata = permute(kdata,[2 3 4 5 1]);
%     kdata = reshape(kdata,[nr*npacc*nz*nttot],nc);
%     [U,S,V] = svd(kdata,'econ');
%     coilinfo = cumsum(diag(S))/sum(diag(S));
%     %     thresh = 1.00;
%     figure(1);clf;
%     plot(1:nc,coilinfo); xlabel('Coils'); ylabel('Normalized Singular Values'); title('PCA Coil Compression');
%     hold on; plot([1 nc],[thresh thresh],'k:','linewidth',1.5);
%     nc = find(coilinfo > thresh,1,'first');
%     fprintf('Number of virtual coils at %.0f%% threshold: %d\n',thresh*100,nc)
%     clear kdata U S
%
%     nVirtualCoils = nc;
%
% end