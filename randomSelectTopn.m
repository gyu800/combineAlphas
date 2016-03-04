function topn = randomSelectTopn( alphas, n )


    [sortedValues,sortIndex] = sort([alphas(:).score],'descend');
    sortedAlphas = alphas(sortIndex);  
    topn = sortedAlphas(1:n);
    
    m = size(alphas,1); 

    if m>n        
        i = 1;
        j = 1;
        maxIndex = zero(n,1);
        while i <=n
            if rand>j/(m+1)
                if isempty(find(maxIndex == sortIndex(j)))
                    maxIndex(i)=sortIndex(j);
                    i=i+1;
                end
            end
            j = j+1;
            if j>m
                j = mod(j,m);                
            end
        end
        topn = alphas(maxIndex); 
    end

end

