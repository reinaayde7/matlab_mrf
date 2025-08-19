clear all, close all, clc
dir = ['E:\POSTDOC_UoM\10_manuscripts\DIP_brainMRF_055T\data_HV_cluster\'];

% %HV 1
% data =  'HV1_ROI_stats_clusterLR.mat';
data =  'HV1_ROI_stats_clusterSLLR.mat';
% data =  'HV1_ROI_stats_clusterDIP.mat';
d{1} = load([dir data]);

% %HV 2
% data =  'HV2_ROI_stats_clusterLR.mat';
data =  'HV2_ROI_stats_clusterSLLR.mat';
% data =  'HV2_ROI_stats_clusterDIP.mat';
d{2} = load([dir data]);

% %HV 3
% data =  'HV3_ROI_stats_clusterLR.mat';
data =  'HV3_ROI_stats_clusterSLLR.mat';
% data =  'HV3_ROI_stats_clusterDIP.mat';
d{3} = load([dir data]);

% %HV 4
% data =  'HV4_ROI_stats_clusterLR.mat';
data =  'HV4_ROI_stats_clusterSLLR.mat';
% data =  'HV4_ROI_stats_clusterDIP.mat';
d{4} = load([dir data]);

%HV 5
% data =  'HV5_ROI_stats_clusterLR.mat';
data =  'HV5_ROI_stats_clusterSLLR.mat';
% data =  'HV5_ROI_stats_clusterDIP.mat';
d{5} = load([dir data]);


%%
T1_WM_LR=[]; T1_WM_LLR=[]; T1_WM_DIP = []; 
T1_GM_LR=[]; T1_GM_LLR=[]; T1_GM_DIP = []; 
T1_CSF_LR=[]; T1_CSF_LLR=[]; T1_CSF_DIP = [];
counter(1) = 1;
counter_gm(1) = 1;
counter_csf(1) = 1;
for v=1:5
    val = d{v}.T1_WM(:);
    val = val(val>20);
    dimension{v}.t1wm_lr = length(val); 
    counter(v+1) = counter(v)+length(val);
    average(v,1) = mean(val);
    stdeviation(v,1) = std(val);
    medianvalue(v,1) = median(val);
    T1_WM_LR = [T1_WM_LR; val];  

    val = d{v}.T1_WM_llr(:);
    val = val(val>20); 
    average(v,2) = mean(val);
    stdeviation(v,2) = std(val);
    medianvalue(v,2) = median(val);
    T1_WM_LLR = [T1_WM_LLR;val];

    val = d{v}.T1_WM_dip(:);
    val = val(val>20); 
    average(v,3) = mean(val);
    stdeviation(v,3) = std(val);
    medianvalue(v,3) = median(val);
    T1_WM_DIP = [T1_WM_DIP;val];


    val = d{v}.T1_GM(:);
    val = val(val>100);
    dimension{v}.t1gm_lr = length(val); 
    counter_gm(v+1) = counter_gm(v)+length(val);
    average(v,4) = mean(val);
    stdeviation(v,4) = std(val);
    medianvalue(v,4) = median(val);
    T1_GM_LR = [T1_GM_LR; val];

    val = d{v}.T1_GM_llr(:);
    val = val(val>100);
    average(v,5) = mean(val);
    stdeviation(v,5) = std(val);
    medianvalue(v,5) = median(val);
    T1_GM_LLR = [T1_GM_LLR;val];

    val = d{v}.T1_GM_dip(:);
    val = val(val>100);
    average(v,6) = mean(val);
    stdeviation(v,6) = std(val);
    medianvalue(v,6) = median(val);
    T1_GM_DIP = [T1_GM_DIP;val];


    val = d{v}.T1_CSF(:);
    val = val(val>100);
    average(v,7) = mean(val);
    stdeviation(v,7) = std(val);
    medianvalue(v,7) = median(val);
    dimension{v}.t1csf_lr = length(val); 
    counter_csf(v+1) = counter_csf(v)+length(val);
    T1_CSF_LR = [T1_CSF_LR; val];

    val = d{v}.T1_CSF_llr(:);
    val = val(val>100);
    average(v,8) = mean(val);
    stdeviation(v,8) = std(val);
    medianvalue(v,8) = median(val);
    T1_CSF_LLR = [T1_CSF_LLR;val];

    val = d{v}.T1_CSF_dip(:);
    val = val(val>100);
    average(v,9) = mean(val);
    stdeviation(v,9) = std(val);
    medianvalue(v,9) = median(val);
    T1_CSF_DIP = [T1_CSF_DIP;val];
   
end
%%
T2_WM_LR=[]; T2_WM_LLR=[]; T2_WM_DIP = []; 
T2_GM_LR=[]; T2_GM_LLR=[]; T2_GM_DIP = []; 
T2_CSF_LR=[]; T2_CSF_LLR=[]; T2_CSF_DIP = [];
counter(1) = 1;
counter_gm(1) = 1;
counter_csf(1) = 1;
for v=1:5
    val = d{v}.T2_WM(:);
    val = val(val>5);
    dimension{v}.t2wm_lr = length(val); 
    counter(v+1) = counter(v)+length(val);
    average(v,1) = mean(val);
    stdeviation(v,1) = std(val);
    medianvalue(v,1) = median(val);
    T2_WM_LR = [T2_WM_LR; val];

    val = d{v}.T2_WM_llr(:);
    val = val(val>5);
    average(v,2) = mean(val);
    stdeviation(v,2) = std(val);
    medianvalue(v,2) = median(val);
    T2_WM_LLR = [T2_WM_LLR;val];

    val = d{v}.T2_WM_dip(:);
    val = val(val>5);
    average(v,3) = mean(val);
    stdeviation(v,3) = std(val);
    medianvalue(v,3) = median(val);
    T2_WM_DIP = [T2_WM_DIP;val];


    val = d{v}.T2_GM(:);
    val = val(val>20);
    dimension{v}.t2gm_lr = length(val); 
    counter(v+1) = counter(v)+length(val);
    average(v,4) = mean(val);
    stdeviation(v,4) = std(val);
    medianvalue(v,4) = median(val);
    T2_GM_LR = [T2_GM_LR; val];

    val = d{v}.T2_GM_llr(:);
    val = val(val>20);
    average(v,5) = mean(val);
    stdeviation(v,5) = std(val);
    medianvalue(v,5) = median(val);
    T2_GM_LLR = [T2_GM_LLR;val];

    val = d{v}.T2_GM_dip(:);
    val = val(val>20);
    average(v,6) = mean(val);
    stdeviation(v,6) = std(val);
    medianvalue(v,6) = median(val);
    T2_GM_DIP = [T2_GM_DIP;val];


    val = d{v}.T2_CSF(:);
    val = val(val>20);
    dimension{v}.t2csf_lr = length(val); 
    counter(v+1) = counter(v)+length(val);
    average(v,7) = mean(val);
    stdeviation(v,7) = std(val);
    medianvalue(v,7) = median(val);
    T2_CSF_LR = [T2_CSF_LR; val];

    val = d{v}.T2_CSF_llr(:);
    val = val(val>20);    
    average(v,8) = mean(val);
    stdeviation(v,8) = std(val);
    medianvalue(v,8) = median(val);
    T2_CSF_LLR = [T2_CSF_LLR;val];

    val = d{v}.T2_CSF_dip(:);
    val = val(val>20);    
    average(v,9) = mean(val);
    stdeviation(v,9) = std(val);
    medianvalue(v,9) = median(val);
    T2_CSF_DIP = [T2_CSF_DIP;val];
end

%%
T1p_WM_LR = T1_WM_LR(T1_WM_LR>20); T1p_WM_LLR = T1_WM_LLR(T1_WM_LLR>20); T1p_WM_DIP = T1_WM_DIP(T1_WM_DIP>20);
T1p_GM_LR = T1_GM_LR(T1_GM_LR>100); T1p_GM_LLR = T1_GM_LLR(T1_GM_LLR>100); T1p_GM_DIP = T1_GM_DIP(T1_GM_DIP>100);
T1p_CSF_LR = T1_CSF_LR(T1_CSF_LR>100); T1p_CSF_LLR = T1_CSF_LLR(T1_CSF_LLR>100); T1p_CSF_DIP = T1_CSF_DIP(T1_CSF_DIP>100);

figure, 
subplot(131), histogram(T1p_WM_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_WM_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_WM_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T1 - WM')
subplot(132), histogram(T1p_GM_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_GM_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_GM_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T1 - WG')
subplot(133), histogram(T1p_CSF_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_CSF_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T1p_CSF_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T1 - CSF')
legend('LR', 'LLR', 'DIP')

T2p_WM_LR = T2_WM_LR(T2_WM_LR>5); T2p_WM_LLR = T2_WM_LLR(T2_WM_LLR>5); T2p_WM_DIP = T2_WM_DIP(T2_WM_DIP>5);
T2p_GM_LR = T2_GM_LR(T2_GM_LR>20); T2p_GM_LLR = T2_GM_LLR(T2_GM_LLR>20); T2p_GM_DIP = T2_GM_DIP(T2_GM_DIP>20);
T2p_CSF_LR = T2_CSF_LR(T2_CSF_LR>20); T2p_CSF_LLR = T2_CSF_LLR(T2_CSF_LLR>20); T2p_CSF_DIP = T2_CSF_DIP(T2_CSF_DIP>20);

figure, 
subplot(131), histogram(T2p_WM_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_WM_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_WM_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T2 - WM')
subplot(132), histogram(T2p_GM_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_GM_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_GM_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T2 - GM')
subplot(133), histogram(T2p_CSF_LR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_CSF_LLR,50, 'EdgeAlpha',0,'Normalization','pdf'), hold on, histogram(T2p_CSF_DIP,50, 'EdgeAlpha',0,'Normalization','pdf'), title('T2 - CSF')
legend('LR', 'LLR', 'DIP')


opt = ['A','B','C','D','E','F','G','H','K','L','M','N','O','P','Q'];

%% T1WM
groups = [repmat("A", dimension{1}.t1wm_lr, 1); 
          repmat("B", dimension{1}.t1wm_lr, 1); 
          repmat("C", dimension{1}.t1wm_lr, 1); 
          repmat("D", dimension{2}.t1wm_lr, 1); 
          repmat("E", dimension{2}.t1wm_lr, 1);
          repmat("F", dimension{2}.t1wm_lr, 1);
          repmat("G", dimension{3}.t1wm_lr, 1); 
          repmat("H", dimension{3}.t1wm_lr, 1);
          repmat("K", dimension{3}.t1wm_lr, 1);
          repmat("L", dimension{4}.t1wm_lr, 1); 
          repmat("M", dimension{4}.t1wm_lr, 1);
          repmat("N", dimension{4}.t1wm_lr, 1);
          repmat("O", dimension{5}.t1wm_lr, 1); 
          repmat("P", dimension{5}.t1wm_lr, 1);
          repmat("Q", dimension{5}.t1wm_lr, 1)];

data = [T1_WM_LR(1:counter(2)-1); T1_WM_LLR(1:counter(2)-1); T1_WM_DIP(1:counter(2)-1); ...
    T1_WM_LR(counter(2):counter(3)-1); T1_WM_LLR(counter(2):counter(3)-1); T1_WM_DIP(counter(2):counter(3)-1);...
    T1_WM_LR(counter(3):counter(4)-1); T1_WM_LLR(counter(3):counter(4)-1); T1_WM_DIP(counter(3):counter(4)-1);...
    T1_WM_LR(counter(4):counter(5)-1); T1_WM_LLR(counter(4):counter(5)-1); T1_WM_DIP(counter(4):counter(5)-1);...
    T1_WM_LR(counter(5):end); T1_WM_LLR(counter(5):end); T1_WM_DIP(counter(5):end);];

x=300:5:700;
labely='T1-WM [ms]';
plotDIST(opt,data,groups,x,labely);

%% T1GM
groups = [repmat("A", dimension{1}.t1gm_lr, 1); 
          repmat("B", dimension{1}.t1gm_lr, 1); 
          repmat("C", dimension{1}.t1gm_lr, 1); 
          repmat("D", dimension{2}.t1gm_lr, 1); 
          repmat("E", dimension{2}.t1gm_lr, 1);
          repmat("F", dimension{2}.t1gm_lr, 1);
          repmat("G", dimension{3}.t1gm_lr, 1); 
          repmat("H", dimension{3}.t1gm_lr, 1);
          repmat("K", dimension{3}.t1gm_lr, 1);
          repmat("L", dimension{4}.t1gm_lr, 1); 
          repmat("M", dimension{4}.t1gm_lr, 1);
          repmat("N", dimension{4}.t1gm_lr, 1);
          repmat("O", dimension{5}.t1gm_lr, 1); 
          repmat("P", dimension{5}.t1gm_lr, 1);
          repmat("Q", dimension{5}.t1gm_lr, 1)];

data = [T1_GM_LR(1:counter(2)-1); T1_GM_LLR(1:counter(2)-1); T1_GM_DIP(1:counter(2)-1); ...
    T1_GM_LR(counter(2):counter(3)-1); T1_GM_LLR(counter(2):counter(3)-1); T1_GM_DIP(counter(2):counter(3)-1);...
    T1_GM_LR(counter(3):counter(4)-1); T1_GM_LLR(counter(3):counter(4)-1); T1_GM_DIP(counter(3):counter(4)-1);...
    T1_GM_LR(counter(4):counter(5)-1); T1_GM_LLR(counter(4):counter(5)-1); T1_GM_DIP(counter(4):counter(5)-1);...
    T1_GM_LR(counter(5):end); T1_GM_LLR(counter(5):end); T1_GM_DIP(counter(5):end);];

x=500:10:1500;
labely='T1-GM [ms]';
plotDIST(opt,data,groups,x,labely);

%% T2WM
groups = [repmat("A", dimension{1}.t2wm_lr, 1); 
          repmat("B", dimension{1}.t2wm_lr, 1); 
          repmat("C", dimension{1}.t2wm_lr, 1); 
          repmat("D", dimension{2}.t2wm_lr, 1); 
          repmat("E", dimension{2}.t2wm_lr, 1);
          repmat("F", dimension{2}.t2wm_lr, 1);
          repmat("G", dimension{3}.t2wm_lr, 1); 
          repmat("H", dimension{3}.t2wm_lr, 1);
          repmat("K", dimension{3}.t2wm_lr, 1);
          repmat("L", dimension{4}.t2wm_lr, 1); 
          repmat("M", dimension{4}.t2wm_lr, 1);
          repmat("N", dimension{4}.t2wm_lr, 1);
          repmat("O", dimension{5}.t2wm_lr, 1); 
          repmat("P", dimension{5}.t2wm_lr, 1);
          repmat("Q", dimension{5}.t2wm_lr, 1)];

data = [T2_WM_LR(1:counter(2)-1); T2_WM_LLR(1:counter(2)-1); T2_WM_DIP(1:counter(2)-1); ...
    T2_WM_LR(counter(2):counter(3)-1); T2_WM_LLR(counter(2):counter(3)-1); T2_WM_DIP(counter(2):counter(3)-1);...
    T2_WM_LR(counter(3):counter(4)-1); T2_WM_LLR(counter(3):counter(4)-1); T2_WM_DIP(counter(3):counter(4)-1);...
    T2_WM_LR(counter(4):counter(5)-1); T2_WM_LLR(counter(4):counter(5)-1); T2_WM_DIP(counter(4):counter(5)-1);...
    T2_WM_LR(counter(5):end); T2_WM_LLR(counter(5):end); T2_WM_DIP(counter(5):end);];

x=0:2:120;
labely='T2-WM [ms]';
plotDIST(opt,data,groups,x,labely);


%% T2GM
groups = [repmat("A", dimension{1}.t2gm_lr, 1); 
          repmat("B", dimension{1}.t2gm_lr, 1); 
          repmat("C", dimension{1}.t2gm_lr, 1); 
          repmat("D", dimension{2}.t2gm_lr, 1); 
          repmat("E", dimension{2}.t2gm_lr, 1);
          repmat("F", dimension{2}.t2gm_lr, 1);
          repmat("G", dimension{3}.t2gm_lr, 1); 
          repmat("H", dimension{3}.t2gm_lr, 1);
          repmat("K", dimension{3}.t2gm_lr, 1);
          repmat("L", dimension{4}.t2gm_lr, 1); 
          repmat("M", dimension{4}.t2gm_lr, 1);
          repmat("N", dimension{4}.t2gm_lr, 1);
          repmat("O", dimension{5}.t2gm_lr, 1); 
          repmat("P", dimension{5}.t2gm_lr, 1);
          repmat("Q", dimension{5}.t2gm_lr, 1)];

data = [T2_GM_LR(1:counter(2)-1); T2_GM_LLR(1:counter(2)-1); T2_GM_DIP(1:counter(2)-1); ...
    T2_GM_LR(counter(2):counter(3)-1); T2_GM_LLR(counter(2):counter(3)-1); T2_GM_DIP(counter(2):counter(3)-1);...
    T2_GM_LR(counter(3):counter(4)-1); T2_GM_LLR(counter(3):counter(4)-1); T2_GM_DIP(counter(3):counter(4)-1);...
    T2_GM_LR(counter(4):counter(5)-1); T2_GM_LLR(counter(4):counter(5)-1); T2_GM_DIP(counter(4):counter(5)-1);...
    T2_GM_LR(counter(5):end); T2_GM_LLR(counter(5):end); T2_GM_DIP(counter(5):end);];

x=0:2:300;
labely='T2-GM [ms]';
plotDIST(opt,data,groups,x,labely);




































function plotDIST(opt,data,groups,x,labely)
    for i=1:length(opt)
         g = strcmp(groups,opt(i));
         vg = data(g);
         pd = fitdist(vg,'kernel');
         dist{i} = pdf(pd,x);
    end
    
   
    fig=figure; 
    fig.Position = [50 50 400 400];
    fig.Color = 'white';
    set(0,"DefaultAxesFontSize",12)
    
    % d1 = [dist{1}',dist{2}'+1e-3,dist{3}'+2e-3];
    tt(1)=0;
    d1 = [dist{1}',dist{2}',dist{3}'];
    M = max(max(d1));
    tt(2)=M;
    plot(d1(:,2),x,'Color','[0 0.4470 0.7410]','LineStyle', '--', 'LineWidth',1.5), hold on
    plot(d1(:,1),x,'Color','[0 0.4470 0.7410]','LineStyle', ':', 'LineWidth',1.5), hold on
    plot(d1(:,3),x,'Color','[0 0.4470 0.7410]','LineStyle', '-', 'LineWidth',1.5), hold on
    
    % d1 = [dist{4}'+M,dist{5}'+1e-3+M,dist{6}'+2e-3+M];
    d1 = [dist{4}'+M,dist{5}'+M,dist{6}'+M];
    M = max(max(d1));tt(3)=M;
    plot(d1(:,2),x,'Color','[0.85 0.3250 0.0980]','LineStyle', '--', 'LineWidth',1.5), hold on
    plot(d1(:,1),x,'Color','[0.85 0.3250 0.0980]','LineStyle', ':', 'LineWidth',1.5), hold on
    plot(d1(:,3),x,'Color','[0.85 0.3250 0.0980]','LineStyle', '-', 'LineWidth',1.5), hold on
    
    % d1 = [dist{7}'+M,dist{8}'+1e-3+M,dist{9}'+2e-3+M];
    d1 = [dist{7}'+M,dist{8}'+M,dist{9}'+M];
    M = max(max(d1));tt(4)=M;
    plot(d1(:,2),x,'Color','[0.929 0.694 0.125]','LineStyle', '--', 'LineWidth',1.5), hold on
    plot(d1(:,1),x,'Color','[0.929 0.694 0.125]','LineStyle', ':', 'LineWidth',1.5), hold on
    plot(d1(:,3),x,'Color','[0.929 0.694 0.125]','LineStyle', '-', 'LineWidth',1.5), hold on
    
    % d1 = [dist{10}'+M,dist{11}'+1e-3+M,dist{12}'+2e-3+M];
    d1 = [dist{10}'+M,dist{11}'+M,dist{12}'+M];
    M = max(max(d1));tt(5)=M;
    plot(d1(:,2),x,'Color','[0.494 0.184 0.556]','LineStyle', '--', 'LineWidth',1.5), hold on
    plot(d1(:,1),x,'Color','[0.494 0.184 0.556]','LineStyle', ':', 'LineWidth',1.5), hold on
    plot(d1(:,3),x,'Color','[0.494 0.184 0.556]','LineStyle', '-', 'LineWidth',1.5), hold on
    
    % d1 = [dist{13}'+M,dist{14}'+1e-3+M,dist{15}'+2e-3+M];
    d1 = [dist{13}'+M,dist{14}'+M,dist{15}'+M];
    M = max(max(d1));tt(6)=M;
    plot(d1(:,2),x,'Color','[0.466 0.674 0.188]','LineStyle', '--', 'LineWidth',1.5), hold on
    plot(d1(:,1),x,'Color','[0.466 0.674 0.188]','LineStyle', ':', 'LineWidth',1.5), hold on
    plot(d1(:,3),x,'Color','[0.466 0.674 0.188]','LineStyle', '-', 'LineWidth',1.5), hold on
    
    grid on
    xticks(tt)
    xlim([0,M])
    xticklabels({'HV1', 'HV2', 'HV3', 'HV4', 'HV5'})
    ylabel(labely,'FontSize', 14, 'FontWeight', 'bold')
end
