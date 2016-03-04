function s_postA = combineAlphas(names,ws,data)
% Alphas
    numAlpha = size(names,1);
    s_postA = 0;
    
    for i = 1:numAlpha
        post = importfile(['20150113_Strategy_Research_Project/All Alpha MAT Files/',char(names(i)), '.postA.mat']);
        postA = post.postA(:,data.ind(1):data.ind(2));
        postA(~isfinite(postA))=0;
        s_postA = s_postA+postA.*ws(i);
    end
end

