% function demo_cpu1_L1PCA
% min_x 0.5 x'*x - ||Ax||_1
clc;clear all;close all;
rand('seed',0); randn('seed',0);
addpath('util','solver','proximal','data');


id_data = [13 23 33];
timeLimits = [2 10 15];

result = [];
for idat = 1:length(id_data)
    G = GetDataSet(id_data(idat),1);
    F = @(x)L1PCA_ComputeTrueObj(x,G);
    for it = 1:5
        fprintf('idata:%d, it:%d\n',idat,it);
        rand('seed',it);
        randn('seed',it);
        x0 = randn(size(G,2),1);
        stopCriteria = 0; % Stop the algorithm when the relative change in F(x) falls below stopCriteria
        timeLimit = 1000;    % Stop the algorithm when the time exceeds timeLimit seconds
        timeLimit = timeLimits(idat);    % Stop the algorithm when the time exceeds timeLimit seconds
        timeIntervel = 0.025;    % Record the objective value in every timeIntervel second
        [x1,his1{it},ts1{it}] = L1PCA_MSCR(G,x0,stopCriteria,timeLimit,timeIntervel);
        [x2,his2{it},ts2{it}] = L1PCA_TolandDual(G,x0,stopCriteria,timeLimit,timeIntervel);
        [x3,his3{it},ts3{it}] = L1PCA_SubGrad(G,x0,stopCriteria,timeLimit,timeIntervel);
        [x4,his4{it},ts4{it}] = L1PCA_CD_SCA(G,x0,stopCriteria,timeLimit,timeIntervel);
        [x5,his5{it},ts5{it}] = L1PCA_CD_SNCA(G,x0,stopCriteria,timeLimit,timeIntervel);
        fprintf('fobj: %f %f %f %f %f, time spent:%f %f %f %f %f \n',F(x1),F(x2),F(x3),F(x4),F(x5), max(ts1{it}), max(ts2{it}), max(ts3{it}), max(ts4{it}), max(ts5{it}) );
    end
    [his1_avg] = GetAvgCell(his1); ts1_avg = GetAvgCell(ts1);
    [his2_avg] = GetAvgCell(his2); ts2_avg = GetAvgCell(ts2);
    [his3_avg] = GetAvgCell(his3); ts3_avg = GetAvgCell(ts3);
    [his4_avg] = GetAvgCell(his4); ts4_avg = GetAvgCell(ts4);
    [his5_avg] = GetAvgCell(his5); ts5_avg = GetAvgCell(ts5);
    
    One = [];
    One.his1 = his1_avg; One.his2 = his2_avg; One.his3 = his3_avg; One.his4 = his4_avg; One.his5 = his5_avg;
    One.ts1 = ts1_avg; One.ts2 = ts2_avg; One.ts3 = ts3_avg; One.ts4 = ts4_avg; One.ts5 = ts5_avg;
    result{idat} = One; save(mfilename,'result')
end



