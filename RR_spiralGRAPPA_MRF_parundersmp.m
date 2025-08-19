function [data_r, ws] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,doCal,cal,ws)

    % data_u: [Nro, TOT_par, Ncoils] -> single time point MRF with all
    % partitions. partition are zeroed when not sampled
    
    % cal: [Nro, par_cal, Ncoils, Nrep] -> pure calibration data. Nrep includes
    % same partition with same spiral projection but difference FA contrast
    
    % pad: ACS Nro block size [-pad RO pad] -> i.e., neighborhood area used to recon RO
    % in the spiral readout dimension 
    
    % kx: kernel RO dimension
    % ky: kernel PAR dimension

    %doCal: runs calibration, can be put on OFF, if ws are available
    %already
    % ws: weights for calibration, if already calculated
    
    
    [nro, npar, ncoils] = size(data_u);
    
    %default settings
    % pad=2;
    % kx = 3;
    % ky=2;
    
    if ky == 2
        p_idx = [-1,1];
    elseif ky == 4
        p_idx = [-3,-1,1,3];
    else %ky=6;
        p_idx = [-5,-3,-1,1,3,5];
    end
    
    %padding to handle recon
    cal_padded = zeros(size(cal,1)+pad*2, size(cal,2), size(cal,3), size(cal,4));
    cal_padded(pad+1:end-pad,:,:,:) = cal;
    data_u_padded = zeros(size(data_u,1)+pad*2, size(data_u,2)+1, size(data_u,3));
    data_u_padded(pad+1:end-pad,1:end-1,:,:) = data_u;

    if doCal == 0 && isempty(ws)
        error("Provide weights for calibration!")
    end
    
    j=1;
    for ro=pad+1:nro %ro: readout point to recon

        if doCal
            cal_i = cal_padded(ro-pad:ro+pad,:,:,:); %calibration [-pad 0 +pad] area around ro point    
    
            % calibration block around RO point
            k=1;
            for r = ceil(kx/2):pad*2+1-floor(kx/2) %indexes in cal_i to be explored with a given kx kernel
                for p =ky:size(cal_i,2)-ky+1 %indexes in cal_i to be recon with a given ky kernel
                    src(:,:,k) = reshape(cal_i(r-floor(kx/2):r+floor(kx/2),p+p_idx,:,:),kx*ky*ncoils, size(cal,4)); %source points
                    trg(:,:,k) = reshape(cal_i(r,p,:,:),1*ncoils,size(cal,4)); %target points
                    k=k+1;
                end
            end
        
            src = reshape(src,kx*ky*ncoils, size(src,2)*size(src,3));
            trg = reshape(trg,1*ncoils, size(trg,2)*size(trg,3));
            ws(:,:,ro) = pinv(src.')*trg.';
        
        
            if j<2
                fprintf('Linear system per each RO: #of Eq. %d ; #of Var. %d ; Ratio of Overdet: %.3f \n',size(src,2), size(src,1), size(src,2) / size(src,1));
            end
            clearvars src trg
        end
        
        % reconstruction: RO point all partitions (NB: accounts only for accel
        % of 2 now
        for p =2:2:npar
            srcr = reshape(data_u_padded(ro-floor(kx/2):ro+floor(kx/2),p+p_idx,:),kx*ky*ncoils,1);
            trgr = (srcr.'*ws(:,:,ro)).'; % target points
            data_u_padded(ro,p,:) = reshape(trgr,[1,ncoils]);
        end
        j=j+1;
    end

    data_r = data_u_padded(pad+1:end-pad,1:end-1,:,:);
end