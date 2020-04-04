% This is the main code of the theoretical limit of distribution system
% identification.

%% We set some hyper parameters
clear;clc;
warning off
caseName = 'case33bw';     % the case name    'case3_dist' 'case33bw'
numSnap = 120;             % the number of snapshot
range.P = 0.6;             % the deviation range of active load 0.6
range.Q = 0.2;             % the deviation range of reactive load to active load 0.2

% the accuracy of measurement device
ratio.P = 0.005;%0.005
ratio.Q = 0.005;
ratio.Vm = 0.0000005;%0.000001; % 0.0000001 0.000001--the maximum error
ratio.Va = 0.005;%0.000005

% if we only compute the bound of admittance matrix
admittanceOnly = false;

% the tolerance of setting a branch to be zero
topoTol = 0.05;

% the enlarge factor to maintain the numerical stability
switch caseName
    case 'case33bw'
        k.G = 5;%1000;
        k.B = 10;%5000;
        k.vm = 100;%100000;
        k.va = 1000;%1000000;
    otherwise
        k.G = 1;
        k.B = 1;
        k.vm = 1;
        k.va = 1;
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
caseDS = caseDS.updateTopo(topoTol, admittanceOnly);
profile off;
% profile viewer;