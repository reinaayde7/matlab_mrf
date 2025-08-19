function [Mx,My,Mz]=dictionary_FISP_withRelaxDecay(r,flip,tr,te,phase,ti,...
    phasetwist,spins,frames,Inv,waitingDuration)


cnt = size(r,1);
nt = length(flip);

if (nt>frames),
    nt = frames;
end;

%df=0;

Mx = zeros(nt*length(te),cnt,'single');
My = zeros(nt*length(te),cnt,'single');
Mz = zeros(nt*length(te),cnt,'single');

phirange = linspace(-phasetwist/2,phasetwist/2,spins);
Z = zeros(3,3,spins);
for s = 1:spins
    Z(:,:,s) = zrot(phirange(s));
end
onemat = ones(1,spins);


for n = 1:cnt
    
    
    T1 = r(n,1)/1000;
    T2 = r(n,2)/1000;
    df = r(n,3);
    
    M0 = repmat([0 0 1]',1,spins);
    M = repmat([0 0 0]',1,spins);
    
    for iEcho = 1:2 %Rudy: run it twice to get to steady state
        %iEcho
        for k = 1:nt
            
            if (k == 1)
                if Inv == 1
                    M = throt(pi,0)*M0;
                end
                [Ati,Bti] = freeprecess(ti,T1,T2,df);
                Bti = Bti*onemat;
                M = Ati*M+Bti;
                
                % spoiling
                for s=1:spins
                    M(:,s)=Z(:,:,s)*M(:,s);
                end
            end;
            
            
            
            
            % alpha pulse
            M = throt(flip(k),phase(k))*M;
            [Ate,Bte] = freeprecess(te,T1,T2,df);
            Bte = Bte*onemat;
            M = Ate*M+Bte;

            if iEcho == 2 %Rudy: run it twice to get to steady state
                Mx(k,n) = mean(M(1,:),2);
                My(k,n) = mean(M(2,:),2);
                Mz(k,n) = mean(M(3,:),2);
            end

            [A,B]=freeprecess(tr(k)-te,T1,T2,df);

            B = B*onemat;
            M=A*M+B;
            for s=1:spins
                M(:,s)=Z(:,:,s)*M(:,s);
            end
            
            
            
        end % timepoints
        
        % relaxation time during the waiting time
        [A,B]=freeprecess(waitingDuration,T1,T2,df);
        B = B*onemat;
        M0=A*M+B;
    end
    %end
    
end % T1 T2 entry
