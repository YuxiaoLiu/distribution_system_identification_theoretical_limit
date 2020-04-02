% This is the main code of the theoretical limit of distribution system
% identification.

%% We set some hyper parameters
clear;
caseName = 'case33bw';     % the case name    'case3_dist' 'case33bw'
numSnap = 120;             % the number of snapshot
range.P = 0.6;             % the deviation range of active load 0.6
range.Q = 0.2;             % the deviation range of reactive load to active load 0.2

% the accuracy of measurement device
ratio.P = 0.05;%0.005
ratio.Q = 0.05;
ratio.Vm = 0.01;%0.000001; % 0.0000001 0.000001--the maximum error
ratio.Va = 0.01;%0.000005

% if we only compute the bound of admittance matrix
admittanceOnly = false;

% the enlarge factor to maintain the numerical stability
switch caseName
    case 'case33bw'
        k.G = 5;%1000;
        k.B = 10;%5000;
        k.vm = 100;%100000;
        k.va = 1000;%1000000;
    case 'case3_dist'
        k.G = 1;
        k.B = 1;
        k.vm = 1;
        k.va = 1;
    otherwise
        k.G = 1;
        k.B = 1;
        k.vm = 50000;
        k.va = 10000000;
end


profile on;
%% We generate the power flow data
caseDS = caseDistributionSystem(caseName, numSnap, range);
caseDS = caseDS.readLoad;
caseDS = caseDS.genOperateData;

%% We evaluate the bound
% set the accuracy of the measurement device, and set whether we have the
% measurement device of a certain state
caseDS = caseDS.setAccuracy(ratio);
caseDS = caseDS.buildFIM(k);
caseDS = caseDS.calBound(admittanceOnly);
% caseDS = caseDS.outputBound;
profile off;
% profile viewer;