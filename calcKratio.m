function k = calcKratio( cumpnl , per)
%Summary of this function goes here
%   The K-ratio is calculated by fitting a linear trend series 
%   to cumulative pnls and estimating the slope and variability of slope
%   Given cumpnl = b * t + a + elpsilon
%   K = b/SE(b) * sqrt(per)/n
%   n is the number of return observations and per corresponds to the expected number of
%  observations in a given calendar year 
%  (e.g. 12 for monthly data, 52 for weekly data, and circa
%  252 for daily data)
%
if nargin < 2
    per = 252;
end

num = size(cumpnl);
x = 1:num(2);
k=zeros(num(1),1);
for i = 1:num(1)
    [slope] = myregr(x,cumpnl(i,:),0);
    k(i) = slope.value/slope.se*sqrt(per)/num(2);
end
k(~isfinite(k )) = 0;

