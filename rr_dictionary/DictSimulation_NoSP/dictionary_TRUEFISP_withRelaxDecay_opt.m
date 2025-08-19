function [Mx,My,Mz]=dictionary_TRUEFISP_withRelaxDecay(r,flip,tr,te,phase,ti,...
    phasetwist,spins,frames,Inv,waitingDuration)

%Rudy: phase twisting is commented out since bSSFP does not add spoilers.
%We commented them out to save sim time.

cnt = size(r,1);
nt = length(flip);

if (nt>frames),
    nt = frames;
end;

%df=0;

Mx = zeros(nt,cnt,'single');
My = zeros(nt,cnt,'single');
Mz = zeros(nt,cnt,'single');

% needed only to dephase after IR pulse
phirange = linspace(-phasetwist/2,phasetwist/2,spins);
Z = zeros(3,3,spins);
for s = 1:spins
    Z(:,:,s) = zrot(phirange(s));
end

onemat = ones(1,spins);


for n = 1:cnt %Rudy debug
    
    
    T1 = r(n,1)/1000; %[s]
    T2 = r(n,2)/1000; %[s]
    df = r(n,3);
    
    M0 = repmat([0 0 1]',1,spins);
    M = repmat([0 0 0]',1,spins);
    
    for iEcho = 1:2 %Rudy: run it twice to get to steady state
        %iEcho
        for k = 1:nt %Rudy debug
            % disp(k)  %Rudy debug
            if (k == 1)
                if Inv == 1
                    M = throt(pi,0)*M0;
                end
                [Ati,Bti] = freeprecess(ti,T1,T2,df);
                Bti = Bti*onemat;
                M = Ati*M+Bti;
                
                %Rudy 2408: spoiler after inversion recovery is kept
                % spoiling
                for s=1:spins
                    M(:,s)=Z(:,:,s)*M(:,s);
                end
            end
            
            
            
            
            % alpha pulse
            M = throt(flip(k),phase(k))*M;
            [Ate,Bte] = freeprecess(te,T1,T2,df);
            Bte = Bte*onemat;
            M = Ate*M+Bte;           

            if iEcho == 2 %Rudy: run it twice to get to steady state
                % M = zrot(-phase(k)+pi/2)*M; %Rudy: rephase signal (ADC phase == RF phase)
                Mx(k,n) = mean(M(1,:),2);
                My(k,n) = mean(M(2,:),2);
                Mz(k,n) = mean(M(3,:),2);
                % M = zrot(phase(k)+pi/2)*M; %Rudy: dephase signal (RF oscillation)
            end
            
            [A,B]=freeprecess(tr(k)-te,T1,T2,df);
            B = B*onemat;
            M=A*M+B;

            % Rudy 240627: spoiler for FISP taken away
            % for s=1:spins
            %     M(:,s)=Z(:,:,s)*M(:,s);
            % end
            
            
            
        end % timepoints
        
        % relaxation time during the waiting time
        [A,B]=freeprecess(waitingDuration,T1,T2,df);
        B = B*onemat;
        M0=A*M+B;
        % disp('done1')
    end
    %end
    % disp('done2')
end % T1 T2 entry
