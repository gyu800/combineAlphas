
clear all; close all;

%data.ws = [10,-0.1,0.1,-10];
%cannot disclose data.mat and alpha data now. 
load('data.mat');
load('alpha1.mat');

[sortedValues,sortIndex] = sort([alpha1(:).score],'descend');
alphas=alpha1(sortIndex);   
%size of subset of alphas for optimizaiton algorithm
n =1000;

%select best nn alphas
topn = selectTopn( alphas, n );
%topn = randomSelectTopn(alpha1', n);  
clear alpha1;
clear alphas;
clear sortedValues;
clear sortIndex;



    
sol=geneticoptimize(topn,data,@costf,1000,20,0.6,0.5, 0.1, 0.01);
sol=geneticoptimize(topn,data,@costf,1000,10,0.6,0.5, 0.1, 0.01);
sol=geneticoptimize(topn,data,@costf,1000,5,0.6,0.5, 0.1, 0.01);
sol=geneticoptimize(topn,data,@costf,1000,2,0.4,0.5,0.2,0.1);
ws = zeros(n,1,'single');
m = 10;
ws(1) = 1;
sol=simplex(topn,ws,data,@costf,0.1,2,1);
sol=simplex(topn,ws,data,@costf,0.1,1000,0.01);

bestAlpha = combineAlphas(dailypnl(topn(sol),:));
bestScore = calcScores(bestAlpha);


bestCumpnl = calcCumPnl(bestAlpha);
bestSharpeR = calcSharpeR(bestAlpha);
bestKratio = calcKratio(bestAlpha);
bestMaxDD= calcMaxDD(bestAlpha);

importfile('C:\Users\Gang\GoogleDrive2\20150113_strategy_data\strategy_data_20150113\data\simvars\dates\20131226_20131231\dates.adj_ammend.mat');    
%t = 1:2013;
%x = datenum(2006, 01, 03) + t;
dates2 = num2str(dates);
x = dates;
for i = 1:size(dates,1)
    x(i) = datenum(str2num(dates2(i,1:4)),str2num(dates2(i,5:6)),str2num(dates2(i,7:8)));
end
    h = figure;
    plot(x(1:1761),bestCumpnl);
    datetick('x','mmm dd yyyy','keeplimits', 'keepticks'); set(gca,'XMinorTick','on','YMinorTick','on')
    ylabel('Cumulative PnL') % y-axis label
    saveas(h,[alphaNames(i,:),'.jpg']);

matfiles = dir(fullfile('C:', 'Users', 'Gang','GoogleDrive2','20150113_Strategy_Research_Project','All Alpha MAT Files','*.mat'));
numAlpha = size(matfiles,1);
for i = 1:numAlpha
    alphaNames(i) = cellstr(matfiles(i).name(1:end-10));
%     loadAlpha3(char(alphaNames(i)));
%     dailypnl(i,:) = pnlSimClose(postA,close);
end
alphaList = alphaNames(topn(sol)')';
