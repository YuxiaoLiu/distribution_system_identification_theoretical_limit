% We generate the boxplot in this figure
addressFigure = '.\plot\case33box.xlsx';

err.gCPS = xlsread(addressFigure, 'gErrCPS');
err.bCPS = xlsread(addressFigure, 'bErrCPS');
bound.g = xlsread(addressFigure, 'gBound');
bound.b = xlsread(addressFigure, 'bBound');
val.gEvalCPS = xlsread(addressFigure, 'gEvalCPS');
val.bEvalCPS = xlsread(addressFigure, 'bEvalCPS');
val.gReal = xlsread(addressFigure, 'gReal');
val.bReal = xlsread(addressFigure, 'bReal');

numBus = size(err.gCPS, 1);
% xlabels = linspace(0,50,numBus);
xlabels = 1:numBus;
boxplot(err.gCPS','width',.5)
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 8);
ax = gca;
ax.XTick = linspace(0,100,numBus);
boxplot(ax, err.gCPS')
xTick = xTick + 1;
hold on;
stem(xTick, bound.g);