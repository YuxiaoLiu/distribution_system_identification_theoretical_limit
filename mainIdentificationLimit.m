% This is the main code of the theoretical limit of distribution system
% identification.

%% We set some hyper parameters
clear;clc;
warning off
caseName = 'case3_dist';     % the case name    'case3_dist' 'case33bw'
numSnap = 20;             % the number of snapshot
range.P = 1.2;               % the deviation range of active load 0.6
range.Q = 0.3;             % the deviation range of reactive load to active load 0.2

% the accuracy of measurement device
ratio.P = 0.001;%0.005
ratio.Q = 0.001;
ratio.Vm = 0.001;%0.0000005;--the maximum error
ratio.Va = 0.001;%0.000005

% if we only compute the bound of admittance matrix
admittanceOnly = false;

% the tolerance of setting a branch to be zero
topoTol = 0.05;

% the prior knowledge of G and B matrix
switch caseName
    case 'case33bw'
        prior.Gmin = 0.3;
        prior.Bmin = 0.3;
    otherwise
        prior.Gmin = 0.1;
        prior.Bmin = 0.1;
end

% % the delta value of FIM matrix
% delta = 0.1;

% the enlarge factor to maintain the numerical stability
switch caseName
    case 'case33bw'
        k.G = 50;%5;%1000;
        k.B = 100;%10;%5000;
        k.vm = 1000;%100;%100000;
        k.va = 10000;%1000;%1000000;
        ratio.Pmin = 0.05; % the minimum P measurement noise compared with the source bus
        ratio.Qmin = 0.05;
    otherwise
        k.G = 1;
        k.B = 1;
        k.vm = 1;
        k.va = 1;
        ratio.Pmin = 0.01; % the minimum P measurement noise compared with the source bus
        ratio.Qmin = 0.01;
end

% profile on;
%% We generate the power flow data
caseDS = caseDistributionSystemMeasure(caseName, numSnap, range);
caseDS = caseDS.readLoad;
caseDS = caseDS.genOperateData;

%% We evaluate the bound
% set the accuracy of the measurement device, and set whether we have the
% measurement device of a certain state
caseDS = caseDS.setAccuracy(ratio);
% profile on
% caseDS = caseDS.buildFIM;
% caseDS = caseDS.calBound;
% caseDS = caseDS.iterateY;
caseDS = caseDS.identifyTopo;
caseDS = caseDS.preEvaluation(prior);
caseDS = caseDS.approximateFIM(k);
caseDS = caseDS.calABound;
% profile off
% profile viewer
% caseDS = caseDS.initValue;
% caseDS = caseDS.identifyMCMCEIO;
% caseDS = caseDS.identifyMCMCEIV;
% caseDS = caseDS.identifyOptNLP;
% caseDS = caseDS.identifyOptGradient;
caseDS = caseDS.identifyOptLMPower;
% caseDS = caseDS.identifyOptNewton;
% caseDS = caseDS.buildFIM(k);
% caseDS = caseDS.updateTopo(topoTol, admittanceOnly);
% profile off;
% profile viewer;