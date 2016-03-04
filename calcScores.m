function alpha = calcScores(dailypnl, ws, corrP)
% scoring system
% factors: fs
% weights: ws
% correlation penalty is optional

    if nargin <3
        alpha.corrP = 0;
    end

    alpha.sharpeR =  calcSharpeR(dailypnl );
    alpha.cumpnl = calcCumPnl(dailypnl);
	alpha.maxDD =  calcMaxDD( alpha.cumpnl );
	alpha.kRatio = calcKratio(alpha.cumpnl,252); 
    

	alpha.score = ws*[alpha.sharpeR;alpha.maxDD;alpha.kRatio;alpha.corrP];

end

