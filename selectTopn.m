function topn = selectTopn( alphas, n )

[sortedValues,sortIndex] = sort([alphas(:).score],'descend');
sortedAlphas = alphas(sortIndex);  
topn = sortedAlphas(1:n);

end

