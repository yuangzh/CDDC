% function demo_cpu6_PhaseRetrival
% min_{s,x} 0.5 ||Gs - xoy||_2^2
% min_{x} - 0.5 x'Ax, s.t. x \in [-1,+1]^n
% clc;clear all;close all;
addpath('util','solver','proximal','data');
rand('seed',0);
randn('seed',0);

id_data = [13 23 33];
timeLimits = [3 3 3];

result = [];
for idat = 1:length(id_data)
    [G,DataStr] = GetDataSet(id_data(idat),2);
    y = GetDataMeasPhaseRetrival(G);
    [Q,normQ] = PhaseRetrivalComputeQ(G,y);
    F = @(x) PhaseRetrival_ComputeTrueObj(x,Q);
    for it = 1:5
        fprintf('idata:%d, it:%d\n',idat,it);
        rand('seed',it);
        randn('seed',it);
        % [s0,~] = eigs(G'*diag(y.*y)*G,1,'LM');[x0] = ComputeOptimalPhaseRetrival(G,s0,y);
        x0 = sign(randn(size(G,1),1));
        stopCriteria = -1; % Stop the algorithm when the relative change in F(x) falls below stopCriteria
        timeLimit = 2;    % Stop the algorithm when the time exceeds timeLimit seconds
        timeLimit = timeLimits(idat);
        timeIntervel = 0.3;    % Record the objective value in every timeIntervel second
        [x1,his1{it},ts1{it}] = PhaseRetrivalAltMin(G,Q,y,x0,stopCriteria,timeLimit,timeIntervel);
        [x2,his2{it},ts2{it}] = PhaseRetrivalGradProj(Q,normQ,x0,stopCriteria,timeLimit,timeIntervel);
        [x3,his3{it},ts3{it}] = PhaseRetrivalGPM(Q,normQ,x0,stopCriteria,timeLimit,timeIntervel);
        [x4,his4{it},ts4{it}] = PhaseRetrivalCD_SCA(Q,normQ,x0,stopCriteria,timeLimit,timeIntervel);
        [x5,his5{it},ts5{it}] = PhaseRetrivalCD_SNCA(Q,normQ,x0,stopCriteria,timeLimit,timeIntervel);
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
