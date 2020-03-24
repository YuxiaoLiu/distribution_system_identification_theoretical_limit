% This is the main code of the theoretical limit of distribution system
% identification.

%% We set some hyper parameters
clc; clear;
caseName = 'case33bw';     % the case name    'case3_dist' 'case33bw'
numSnap = 30;              % the number of snapshot
range.P = 0.6;             % the deviation range of active load
range.Q = 0.2;             % the deviation range of reactive load to active load

profile on;
%% We generate the power flow data
caseDS = caseDistributionSystem(caseName, numSnap, range);
caseDS = caseDS.readLoad;
caseDS = caseDS.genOperateData;

%% We evaluate the bound
% set the accuracy of the measurement device, and set whether we have the
% measurement device of a certain state
caseDS = caseDS.setAccuracy;
caseDS = caseDS.buildFIM;
profile off;
profile viewer;