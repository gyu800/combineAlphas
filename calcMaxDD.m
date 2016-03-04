
function maxDD = calcMaxDD(cumpnl)
    %cumpnl(numAlphas,days);

    maxDD = zeros(size(cumpnl,1),1);
    for x = 1:size(cumpnl,1)
    
        n = realmax;
        high = 0;
        low = n;
        for i = 1:size(cumpnl,2)
            if cumpnl(x,i) > high
                high = cumpnl(x,i);
                low = n;
            end

            if cumpnl(x,i) < low
                low = cumpnl(x,i);
                if low < high
                    dd = (high - low)/abs(high);
                    if dd > maxDD(x)
                        maxDD(x) =dd;
                    end
                end
            end
        end     
    end
    
    %maxDD = maxDD./max(cumpnl(:,end))*100;
    %DD is percentage of total pnl.
    maxDD (~isfinite(maxDD )) = 0;
    %maxDD (maxDD>1) = 1;
end