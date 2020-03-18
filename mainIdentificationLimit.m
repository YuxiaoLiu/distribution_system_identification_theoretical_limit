% This is the main code of the theoretical limit of distribution system
% identification.
%% We first set some hyper parameters
clc; clear;
caseName = 'case33bw';     % the case name
numSnap = 50;              % the number of snapshot
%% Then we generalize the power flow data
caseDS = caseDistributionSystem(caseName, numSnap);