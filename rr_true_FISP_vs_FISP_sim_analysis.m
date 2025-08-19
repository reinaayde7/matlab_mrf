clear all, close all, clc
%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
fatxtfile = [folderTXT '\FISP_FA_Body.txt'];
% fatxtfile = [folderTXT '\FISP_FA_body_x1_5.txt'];
% fatxtfile = [folderTXT '\FISP_FA_body_x2.txt'];
phtxtfile = [folderTXT '\FISP_PH.txt'];


log.FA = importdata(fatxtfile);
% log.tr0 = log.rawinfo.TR(1)/1000/1000;
log.TR = 14;%importdata(trtxtfile)*1e-6 + log.tr0;
log.TE = 1.8;%importdata(tetxtfile);
% log.PH = importdata(phtxtfile);

% plot the parameters
figure('name','Acquisition Parameters');
plot(log.FA(1:600),'LineWidth',2);
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');
%%
FA = log.FA(1:600);
te(1) = 1.8;
tr = ones(600,1)*14;
for i = 2:600
    te(i) = ( tr(i-1)-te(i-1) )*sin( FA(i-1)/2) / sin(FA(i)/2 );
end

%%
figure, plot(te)

%%


%% Dictionary handler: 

% dictionary_name = 'dict_trueFISP_TE1_8_TR13_FAbody15_Nex600_b0100';
dictionary_name = 'dict_trueFISP_TE1_8_TR13_FAbody10_Nex600_b0100';
dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
dict_truefisp = dfile.dict; r_truefisp = dfile.r; dict_log_truefisp = dfile.log; 
clearvars dfile        

%%
% dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody15_Nex600';
dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody10_Nex600';
dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
dict_fisp = dfile.dict; r_fisp = dfile.r; dict_log_fisp = dfile.log;
clearvars dfile

%%
dict_fisp = squeeze(dict_fisp);
dict_truefisp_b0 = dict_truefisp(:,r_truefisp(:,3)==0);

%%
t1 = 490;
t2=50;

b0 = 0;
r_idx_truefisp0 = find(r_truefisp(:,1) == t1 & r_truefisp(:,2) == t2 & r_truefisp(:,3) == b0);
b0 = 20;
r_idx_truefisp20 = find(r_truefisp(:,1) == t1 & r_truefisp(:,2) == t2 & r_truefisp(:,3) == b0);
b0 = 40;
r_idx_truefisp40 = find(r_truefisp(:,1) == t1 & r_truefisp(:,2) == t2 & r_truefisp(:,3) == b0);

r_idx_fisp = find(r_fisp(:,1) == t1 & r_fisp(:,2) == t2);

%%

s_fisp = dict_fisp(:,r_idx_fisp);
powerfisp = mean(abs(s_fisp).^2);
s_truefisp0 = dict_truefisp(:,r_idx_truefisp0);
powertruefisp0 = mean(abs(s_truefisp0).^2);
s_truefisp20 = dict_truefisp(:,r_idx_truefisp20);
powertruefisp20 = mean(abs(s_truefisp20).^2);
s_truefisp40 = dict_truefisp(:,r_idx_truefisp40);
powertruefisp40 = mean(abs(s_truefisp40).^2);

[min20,max20,power20] = calc_ranges(t1,t2,r_truefisp, dict_truefisp,[-20:1:20]);
[min40,max40,power40] = calc_ranges(t1,t2,r_truefisp, dict_truefisp,[-40:1:40]);
%%


x = 1:600;
figure,
set(gcf, 'Color', 'white')
fill([x fliplr(x)], [min40; flipud(max40)], [0 0 1], 'EdgeColor', 'none','FaceAlpha',0.1), hold on
fill([x fliplr(x)], [min20; flipud(max20)], [0 0 1], 'EdgeColor', 'none','FaceAlpha',0.3)
hold on, plot(x,abs(s_truefisp0), 'b', 'LineWidth',2), hold on
hold on, plot(x,abs(s_truefisp20), 'b--', 'LineWidth',2), hold on
hold on, plot(x,abs(s_truefisp40), 'b-.', 'LineWidth',2), hold on
plot(abs(s_fisp), 'r', 'LineWidth',2)
a = sprintf(['pSSFP B0: [-40,+40] Hz, avg.pow: %.2f' char(177) '%.2f'], mean(power40)*1e2, std(power40)*1e2);
b = sprintf(['pSSFP B0: [-20,+20] Hz, avg.pow: %.2f' char(177) '%.2f'], mean(power20)*1e2, std(power20)*1e2);
c = sprintf('pSSFP B0: 0 Hz, avg.pow: %.2f', mean(powertruefisp0)*1e2);
d = sprintf('pSSFP B0: 20 Hz, avg.pow: %.2f', mean(powertruefisp20)*1e2);
e = sprintf('pSSFP B0: 40 Hz, avg.pow: %.2f', mean(powertruefisp40)*1e2);
f = sprintf('FISP avg.pow: %.2f', mean(powerfisp)*1e2);
legend(a,b,c,d,e,f)
title(sprintf('T1: %.0f ms, T2: %.0f ms', t1, t2))
ylim([0, 0.45])
xlabel('Image Time Points')
ylabel('|S|')


%% calculation of power ratio as function of T1 and T2, as well as B0
B = 35;
idx = r_truefisp(:,3)==B;
dict_truefisp_B = dict_truefisp(:,idx);
for j=1:size(dict_fisp,2)
    p_fisp(j) = mean(abs( dict_fisp(:,j)).^2);
    p_truefisp(j) = mean(abs( dict_truefisp_B(:,j)).^2);
end

%%
rp = p_truefisp./p_fisp;
%%
x = r_fisp(:,1); y = r_fisp(:,2);
z = rp;

i = find(z<1.02 & z>0.98);


r1 = r_truefisp(idx,:); d1 = dict_truefisp(:,idx); %B0
jjj =  r1(:,1) < 600 & r1(:,1) > 300; 
r2 = r1(jjj,:); d2 = d1(:,jjj); %T1
jjjj = r2(:,2) < 100;
r3 = r2(jjjj,:); d3 = d2(:,jjjj); %T2
p_wm_truefisp = mean(abs(d3).^2);

jjj =  r_fisp(:,1) < 600 & r_fisp(:,1) > 300; 
r1 = r_fisp(jjj,:); d1 = dict_fisp(:,jjj); %T1
jjjj = r1(:,2) < 100;
r2 = r1(jjjj,:); d2 = d1(:,jjjj); %T2
p_wm_fisp = mean(abs(d2).^2);

rp_wm = mean(p_wm_truefisp./p_wm_fisp);

r1 = r_truefisp(idx,:); d1 = dict_truefisp(:,idx); %B0
jjj =  r1(:,1) < 1000 & r1(:,1) > 600; 
r2 = r1(jjj,:); d2 = d1(:,jjj); %T1
jjjj = r2(:,2) < 200 & r2(:,2) > 50; 
r3 = r2(jjjj,:); d3 = d2(:,jjjj); %T2
p_gm_truefisp = mean(abs(d3).^2);

jjj =  r_fisp(:,1) < 1000 & r_fisp(:,1) > 600; 
r1 = r_fisp(jjj,:); d1 = dict_fisp(:,jjj); %T1
jjjj = r1(:,2) < 200 & r1(:,2) > 50; 
r2 = r1(jjjj,:); d2 = d1(:,jjjj); %T2
p_gm_fisp = mean(abs(d2).^2);

rp_gm = mean(p_gm_truefisp./p_gm_fisp);
%%
figure, scatter(x,y,100, z, 'filled'), hold on

xc = x(i);
yc = y(i);
[xc_s, jj] = sort(xc);
plot(xc_s, yc(jj),'k', 'LineWidth',3)
colormap('jet'), colorbar
xlabel('T1 [ms]'), ylabel('T2 [ms]')
title(['pSSFP/FISP power ration - B0: ' num2str(B) 'Hz'])
xlim([0 1500])
ylim([0 250])
caxis([0 3])
rectangle('Position', [300 0 300 100], 'EdgeColor', 'w', 'LineWidth', 2);
text(400, 50, num2str(round(rp_wm,2)), 'FontSize',12, 'Color', 'w')
rectangle('Position', [600 50 400 150], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
text(750, 125, num2str(round(rp_gm,2)), 'FontSize',12, 'Color', [0.5 0.5 0.5])






%%
function [minV,maxV,powerV] = calc_ranges(t1,t2,r_truefisp, dict_truefisp,b0_ranges)
    i=1;
    for n = 1:numel(b0_ranges)
        b0 = b0_ranges(n);
        r_idx_truefisp = find(r_truefisp(:,1) == t1 & r_truefisp(:,2) == t2 & r_truefisp(:,3) == b0);
        s_truefisp_b0var(:,i) = dict_truefisp(:,r_idx_truefisp);
        i=i+1;
    end
    minV = min(abs(s_truefisp_b0var),[],2);
    maxV = max(abs(s_truefisp_b0var),[],2);
    powerV = mean(abs(s_truefisp_b0var).^2);
    % figure, hist(powerV*1e2,30)
end
















