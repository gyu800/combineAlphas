
function sharpeR = calcSharpeR(dailypnl)

    dailypnl(~isfinite(dailypnl)) = 0;
    %dailypnl(numAlphas,days);
    
    % dolar neutral so no substraction of risk-free return.
    sharpeR = mean(dailypnl,2)./std(dailypnl,0,2);
    sharpeR(~isfinite(sharpeR)) = 0;
end