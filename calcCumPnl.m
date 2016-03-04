
function cumpnl = calcCumPnl(dailypnl)
    %dailypnl(numAlphas,days);
    
    cumpnl = dailypnl;
    for n = 2:size(cumpnl,2)
        cumpnl(:,n) = cumpnl(:,n) + cumpnl(:,n-1);
    end
    
end