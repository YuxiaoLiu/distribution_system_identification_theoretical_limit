% This is the main code of the theoretical limit of distribution system
% identification.
%% We first set some hyper parameters
clc; clear;
caseName = 'case33bw';     % the case name
numSnap = 50;              % the number of snapshot
range.P = 0.6;             % the deviation range of active load
range.Q = 0.2;             % the deviation range of reactive load to active load
%% Then we generalize the power flow data
caseDS = caseDistributionSystem(caseName, numSnap, range);