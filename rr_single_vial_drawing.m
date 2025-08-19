% investigate vial #6 T1 790, t2 216

meas = avg_fp(6,:);
[t1e,t2e,b0e, iP] = patternmatchDEBUG(meas,1,r,0,dict,dictnorm,1,1);

%%
possible_sol = find(abs(iP)>0.6337);
figure, plot(abs(iP(possible_sol))), ylabel('Inner Prod'), xlabel('possible solutions')

figure, 
subplot(311), plot(r(possible_sol,1)), ylabel('T1 [ms]')
subplot(312), plot(r(possible_sol,2)), ylabel('T2 [ms]')
subplot(313), plot(r(possible_sol,3)), ylabel('B0 [Hz]'), xlabel('possible solutions')

%%
t1lookup= [790, 760];
t2lookup= [216, 135];
b0lookup= [10,10];

[idxsLU, valuesLU, fpLU] = dict_entry_lookup(dict, r, t1lookup, t2lookup, b0lookup);


%%
xx = imag(meas) + j*real(meas);
xxNorm = sqrt(sum(xx.*conj(xx),2));

figure, 
subplot(311)
plot(abs(xx/xxNorm),'b'), hold on
plot(abs(squeeze(fpLU(:,1)))/dictnorm(idxsLU(1)), 'LineWidth',2), hold on
plot(abs(squeeze(fpLU(:,2)))/dictnorm(idxsLU(2)), 'LineWidth',2),
ylabel('abs')
subplot(312)
plot(real(xx/xxNorm),'b'), hold on
plot(real(squeeze(fpLU(:,1)))/dictnorm(idxsLU(1)), 'LineWidth',2), hold on
plot(real(squeeze(fpLU(:,2)))/dictnorm(idxsLU(2)), 'LineWidth',2),
ylabel('real')
subplot(313)
plot(imag(xx/xxNorm),'b'), hold on
plot(imag(squeeze(fpLU(:,1)))/dictnorm(idxsLU(1)), 'LineWidth',2), hold on
plot(imag(squeeze(fpLU(:,2)))/dictnorm(idxsLU(2)), 'LineWidth',2),
xlabel('Nex'), ylabel('imag')


%%

v1 = squeeze(fpLU(:,1))/dictnorm(idxsLU(1));
v2 = squeeze(fpLU(:,2))/dictnorm(idxsLU(2));
corr = v1'*v2;

abs(corr)


%%
t1lookup= [230];
t2lookup= [50];
b0lookup= [0];

[idxsLU, valuesLU, fpLU] = dict_entry_lookup(dict, r, t1lookup, t2lookup, b0lookup);

signal = dict(:,idxsLU);     % conversion needed to use old code based on invivo data;
[tt1,tt2,bb0,iP] = patternmatchDEBUG(signal',1,r,0,dict,dictnorm,1,1);

%%
figure, plot(abs(iP))