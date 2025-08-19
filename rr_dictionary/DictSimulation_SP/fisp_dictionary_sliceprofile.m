function [dict] = fisp_dictionary_sliceprofile(r,flip,tr,ti,x,grad,rf,F,t,ResidualGrad)
cnt = size(r,1);
NumOfFrames = length(flip);
dict = complex(zeros(NumOfFrames,cnt,'single'),zeros(NumOfFrames,cnt,'single'));

for jj = 1:cnt,
    
    if mod(jj,500)==0
        fprintf('MRF Dictionary: %.0f / %.0f\n',lEntry,cnt);
    end
    
    t1 = r(jj,1);
    t2 = r(jj,2);
    NumOfSpins = length(x);
    mx = zeros(NumOfSpins,1);
    my = zeros(NumOfSpins,1);
    mz = -1*ones(NumOfSpins,1);
    
    
    [mx,my,mz] = bloch(0,0,ti/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
    for ii = 1:NumOfFrames,
        
        theta1 = flip(ii)/180*pi;
        
        A1 = (theta1./((42.58*100)*2*pi*F*((t(2)-t(1))/(1000*1000))));%T->gauss
        rf_scaled = (A1*rf);
        
        
        [mx,my,mz] = bloch(rf_scaled,grad/10,(t(2)-t(1))/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        
        [mx,my,mz] = bloch(0,0,50/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        msig=mx+1i*my;
        [mx,my,mz] = bloch(0,0,(tr(ii) - t(end)-50)/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        %additional dephasing gradient
        [mx,my,mz] = bloch(0,ResidualGrad(ii)/10,10/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        
        
        dict(ii,jj) = sum(msig);
        
        
    end
end
