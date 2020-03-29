% This is the main code of the theoretical limit of distribution system
% identification.

%% We set some hyper parameters
clc; clear;
caseName = 'case33bw';     % the case name    'case3_dist' 'case33bw'
numSnap = 120;             % the number of snapshot
range.P = 0.6;             % the deviation range of active load 0.6
range.Q = 0.2;             % the deviation range of reactive load to active load 0.2

% the accuracy of measurement device
ratio.P = 0.005;%0.005
ratio.Q = 0.005;
ratio.Vm = 0.0001; % 0.0000001 0.000001--the maximum error
ratio.Va = 0.0001;

% if we only compute the bound of admittance matrix
admittanceOnly = true;

profile on;
%% We generate the power flow data
caseDS = caseDistributionSystem(caseName, numSnap, range);
caseDS = caseDS.readLoad;
caseDS = caseDS.genOperateData;

%% We evaluate the bound
% set the accuracy of the measurement device, and set whether we have the
% measurement device of a certain state
caseDS = caseDS.setAccuracy(ratio);
caseDS = caseDS.buildFIM;
caseDS = caseDS.calBound(admittanceOnly);
% caseDS = caseDS.outputBound;
profile off;
profile viewer;