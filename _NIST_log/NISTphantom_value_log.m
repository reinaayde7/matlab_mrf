% create .mat file with NIST values

% values from Jesus's table sent by email on 240627 as a png 
NIST.T1array.T1.freemax = [1716.6, 1269.2, 950.0, 702.2, 509.7, 365.9, 260.5, 191.5, 133.0, 94.2, 66.9, 47.3, 34.0, 24.3];
NIST.T1array.T2.freemax = [1500.0, 1234.3, 946.1, 683.0, 501.1, 361.5, 257.1, 187.9, 131.3, 92.9, 65.7, 47.3, 33.5, 24.1];

% values from Jesus's table sent by email on 240627 as a png 
% NIST.T2array.T1.freemax = [1985.9, 1761.5, 1525.0, 1246.9, 1037.0, 803.0, 611.2, 468.3, 348.5, 233.4, 181.8, 127.2, 94.5, 66.7];
% NIST.T2array.T2.freemax = [1043.7, 789.3, 580.8, 416.9, 306.7, 216.1, 152.1, 111.3, 78.5, 52.1, 39.0, 26.9, 19.9, 14.2];

% values measured Jun 2024 by Tom
NIST.T2array.T1.freemax = [1990, 1756, 1519, 1241, 1018, 790.4, 600.2, 458.7, 345.9, 230.6, 179.4, 125.7, 92.8, 65.2];
NIST.T2array.T2.freemax = [1066, 797.5, 579.8, 415.5, 307.6, 215.9, 149, 108.4, 78.28, 49.8, 38.1, 25.8, 19.1, 13];

savedir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
filename = ['NIST_log' date];
save([savedir filename],'NIST','-v7.3');