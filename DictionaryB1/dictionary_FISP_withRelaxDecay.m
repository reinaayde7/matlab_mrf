function [Mx,My,Mz]=dictionary_FISP_withRelaxDecay(r,flip,tr,te,phase,ti,...
    phasetwist,spins,frames,Inv,waitingDuration)


cnt = size(r,1);
nt = length(flip);

if (nt>frames),
    nt = frames;
end;

df=0;

Mx = zeros(nt,cnt,'single');%Originally zeros(nt*length(te),cnt,'single');, Jesus changed it
My = zeros(nt,cnt,'single');
Mz = zeros(nt,cnt,'single');


phirange = linspace(-phasetwist/2,phasetwist/2,spins);
Z = zeros(3,3,spins);
for s = 1:spins
    Z(:,:,s) = zrot(phirange(s));
end
onemat = ones(1,spins);


for n = 1:cnt
       
    T1 = r(n,1)/1000;
    T2 = r(n,2)/1000;
    b1map = r(n,3);
    
    M0 = repmat([0 0 1]',1,spins,length(te)+1); 
    M = repmat([0 0 0]',1,spins);
    
    for iEcho = 1:length(te)
        %iEcho
        for k = 1:nt
            
            if (k == 1)
                if Inv == 1
                    M = throt(pi,0)*M0(:,:,iEcho);
                end
                [Ati,Bti] = freeprecess(ti,T1,T2,0);
                Bti = Bti*onemat;
                M = Ati*M+Bti;
                
                % spoiling
                for s=1:spins
                    M(:,s)=Z(:,:,s)*M(:,s);
                end
            end;
            
            % alpha pulse
            %disp('flip(k): ')
            %disp(num2str(flip(k)))
            %disp('B1map :')
            %disp(num2str(((b1map))))
            %disp('flip(k)*B1map :')
            %disp(num2str((flip(k)*(b1map))))
            M = throt(flip(k)*(b1map),phase(k))*M; %Jesus added *b1map
            [Ate,Bte] = freeprecess(te(iEcho),T1,T2,0);
            Bte = Bte*onemat;
            M = Ate*M+Bte;
            Mx(nt*(iEcho-1)+k,n) = mean(M(1,:),2);
            My(nt*(iEcho-1)+k,n) = mean(M(2,:),2);
            Mz(nt*(iEcho-1)+k,n) = mean(M(3,:),2);
            [A,B]=freeprecess(tr(k)-te(iEcho),T1,T2,0);
            B = B*onemat;
            M=A*M+B;
            for s=1:spins
                M(:,s)=Z(:,:,s)*M(:,s);
            end
            
            
            
        end % timepoints
        
        % relaxation time during the waiting time
        [A,B]=freeprecess(waitingDuration,T1,T2,0);
        B = B*onemat;
        M0(:,:,iEcho+1)=A*M+B;
    end
    %end
    
end % T1 T2 entry
