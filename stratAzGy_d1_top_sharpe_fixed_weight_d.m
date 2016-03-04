% Strategy template for Strategy Trainee Candidates.
% Copyright 2015 Trexquant Investment LP - All Rights Reserved.
% Redistribution without written consent from Trexquant is prohibited.

% TODO: Put your initials in place of "Yi" below.
%       For example if your name is John Doe the file would be "stratAzJd_d1...".
%       Then change the name of this file itself to match the function name,
%       i.e., resave the file as "stratAzJd_d1...m".
function [preA postA pre_simres post_simres alpha_weights_matrix] = ...
			stratAzGy_d1_top_sharpe_fixed_weight(data, alpha_cube, varargin)

params = inputParser;

% TODO: tighten constraints below if desired.  Never loosen them.  Tightening constraints may hurt performance.
params.addParamValue('pos_limit', .005); % the maximum absolute position for a single stock in the portfolio (<= 0.005)
params.addParamValue('max_trd', 0.01); % may only trade 1% of vol20 per day (<= 0.01)
params.addParamValue('max_pos', 0.30); % may only hold a maximum absolute position of 30% of vol20*price (<= 0.30)

% Do not edit the parameters below.
params.addParamValue('booksize', 500e6); % the long-side size of the portfolio
params.addParamValue('universe', 'TOP2000'); % stock universe to trade after applying all post processing
params.addParamValue('stockpnl', 0, @(x)(x==0)||(x==1)); % saves pnl per stock (high disk usage if enabled)
params.addParamValue('should_position_pnl', 0, @(x)(x==0)||(x==1)); % saves position pnl per stock (high disk usage if enabled)
params.addParamValue('should_execution_pnl', 0, @(x)(x==0)||(x==1)); % saves execution pnl per stock (high disk usage if enabled)
params.addParamValue('should_slippage_pnl', 0, @(x)(x==0)||(x==1)); % saves slippage pnl per stock (high disk usage if enabled) 
params.addParamValue('should_dividend_pnl', 0, @(x)(x==0)||(x==1)); % saves dividend pnl per stock (high disk usage if enabled)
params.addParamValue('delay', 1,  @(x)(x==0)||(x==1)||(x==2)); % target is ready in the morning (no intraday data)
params.addParamValue('locate_utilization_rate', .1); % the percent of the locate availability we may actually use
params.addParamValue('max_rebate_rate', -5); % the maximum (negative) rebate rate we are willing to pay to short, in basis points
params.addParamValue('startdate', 20100101); % start date for simulation chart in report email
params.addParamValue('yesterday_preA', []); % production-related parameter
params.addParamValue('meta', 0);
params.parse(varargin{:});
pos_limit = params.Results.pos_limit;
max_trd = params.Results.max_trd;
max_pos = params.Results.max_pos;
booksize = params.Results.booksize;
universe = params.Results.universe;
stockpnl = params.Results.stockpnl;
should_position_pnl = params.Results.should_position_pnl;
should_execution_pnl = params.Results.should_execution_pnl;
should_slippage_pnl = params.Results.should_slippage_pnl;
should_dividend_pnl = params.Results.should_dividend_pnl;
delay = params.Results.delay;
locate_utilization_rate = params.Results.locate_utilization_rate;
max_rebate_rate = params.Results.max_rebate_rate;
startdate = params.Results.startdate;
yesterday_preA = params.Results.yesterday_preA;
meta = params.Results.meta;

meta_data = struct();

meta_data.alpha_name = 'stratAzGy_d1_top_sharpe_fixed_weight'; % TODO: alpha_name must be identical to filename and strategy function name.
meta_data.authors = {'Adam Zelinski'}; %%%%% 'Gary Yu'}; % TODO: fill in your name.  Adam's information is listed first so he receives error emails and copies of simulation reports.
meta_data.username = {'azelinski'}; %%%%%, 'gyu'}; % TODO: fill in your abbreviated name
meta_data.strat_id = [];
meta_data.author_emails = {'adam.zelinski@trexquant.com'}; %%%%%% ,'yu.gary1@gmail.com'}; % , 'Your.Email@Your.Email.Host.com'}; % TODO: fill in your e-mail.
meta_data.region_developed_for = 'USA';
meta_data.expected_running_time = 60;
meta_data.birthday = 20131231; % TODO: if you use all data <= YYYYMMDD to create a list of fixed alpha weights, then "birthday" equals YYYYMMDD.
                               % See detailed strategy description below to understand concept of birthday.
meta_data.keywords = {'delay1', 'rw'}; % TODO: put keywords here, e.g., rw = rolling window
meta_data.description = ['Simplex look back 500 days']; ...]; % TODO: add short description
meta_data.use_shared_memory_alpha_cube = 1;
meta_data.variables_used = {'ret1','vol20','locates_ms_avshs','locates_ms_rebate','industry'}; % TODO: many variables are already in the data struct.  If you need others, simply list them here.
                                     % TODO: any variable in the simvars directory you downloaded is loadable by specifying it in the cell array.
meta_data.test_startdates = [20060302 20100101];

if(meta)
	preA = meta_data;
	disp('not running -- returning meta_data');
	return;
end

% TODO: understand the two key structures explained below.
%
% (1) data: this is a structure of any data in the simvars folder.  
%       By default it will have close, volume, etc. as fields.  
%       Any simvar listed in meta_data.variables_used will also be in the data structure.
%       All data is adjusted for relevant corporate actions (e.g., stock splits).
%
%           data.close: 2-D matrix, numstocks x numdays, closing price per stock per day
%           data.volume: 2-D matrix, numstocks x numdays, shares traded per stock per day
%
% (2) alpha_cube: this is a read-only structure containing all information for all 4,000+ alphas.
%                 attempting to write to any field of this data structure is prohibited and will result in the strategy erroring out.
%
%         simres: {1x4932 cell} -- there is a simres for each alpha
%          names: {1x4932 cell} -- there is a name for each alpha
%         values: [4932x5777x1236 single] -- the 2-D numstocks x numdays matrices for each alpha are stacked up into a 3-D matrix
%
%     * alpha_cube.simres{1}:
%			dates: [1236x1 int32]
%		 longsize: [1236x1 single]
%		shortsize: [1236x1 single]
%		  longstk: [1236x1 single]
%		 shortstk: [1236x1 single]
%		 dailypnl: [1236x1 single]
%			   ir: [1236x1 single]
%			  ret: [1236x1 single]
%			  tvr: [1236x1 single]
%		 dailytvr: [1236x1 single]
%		centpersh: [1236x1 single]
%			   dd: [1236x1 single]
%			delay: 1
%		  liqsize: [1236x1 single]
%
%     * alpha_cube.names{1}:
%
%          ALPHA_1389
%
%     * squeeze(alpha_cube.values(1, :, :)) is the full 2-D matrix for ALPHA_1389.
%

if(length(universe) > 0)
	data.valids = data.(lower(universe));
end

numdates = numel(data.dates);
numalphas = numel(alpha_cube.names);
numstocks = data.numstocks;

% TODO: each day, your goal is to designate a weight per alpha.
alpha_weights_matrix = zeros(numalphas, data.numdates, 'single');

disp('size(alpha_weights_matrix):') 
disp(num2str(size(alpha_weights_matrix)))

% TODO: read and make sure you understand the fixed-weight strategy implemented below.

%%%%%%%%%%%% Simple Fixed-Weight Strategy %%%%%%%%%%%%

% Put all simres PnL in 1 place.
alpha_PnL = [];
for ai = 1 : numalphas,
	alpha_PnL = [alpha_PnL; alpha_cube.simres{ai}.dailypnl.'];
end
disp(['size(alpha_PnL, 1): ' num2str(size(alpha_PnL, 1))])
disp(['size(alpha_PnL, 2): ' num2str(size(alpha_PnL, 2))])

% TODO: Understand the fixed weight strategy example given below.
% By "fixed weight", we mean we will fill in the columns of the alpha_weights_matrix 
% (the weights to apply to the alphas on a given day) with a single set of weights.
% The strategy will use any data structure and alpha_cube structure information up to 20131231.
% It will generate a list of weights to apply to each alpha.
% Then it will apply these weights all days in the future.
% Fixed-weight strategies have very clear in-sample vs. out-of-sample periods:
%   In-sample: all dates < 20140102
%   Out-of-sample: all dates >= 20140102

di_start = find(data.dates == 20131231);
back_day = find(data.dates == data.dates(max(1, di_start - 500))); % prevent user from looking back before data exists
num_top_alpha = 50;

% On di_start, we look backward by up to 500 days and compute each alpha's Sharpe ratio.
% The top 50 alphas will receive positive weight.  All other ~4900 alphas will receive zero weight.
% The weights will be proportional to the Sharpe ratio of each of the 50 alphas.
% The same weights will be applied for all days before and after di_start.
%   * Backwards in time (in-sample), to judge how the strategy works in-sample.
%   * Forwards in time (out of sample), to judge how the strategy works out of sample.
% Generally, weights applied in the in-sample period have good PnL per day.
% It is more challenging to generate strong PnL during the out-of-sample period.
%
% *** To be considered for live trading, a strategy should have a good Sharpe ratio and 
%     return in the out-of-sample period.  That is your goal. ***
%
% *** Never "look ahead" at future data, or the result is forward biased, and not tradeable.
%     In other words, on day index di, never use data.close(:, di+1) or alpha_cube.values(:, :, di+1).
%     On day index di, only use data and alpha_cube information up to and including day index di. ***
%

for di = di_start,
	% Display progress.
	disp(['di = ' num2str(di)])
	
	% Looking backward from day index di, compute Sharpe of each alpha.
	PnL_mean_train = mean(alpha_PnL(:, (di - back_day + 1):di), 2);
	PnL_var_train = var(alpha_PnL(:, (di - back_day + 1):di), 0, 2);

	PnL_sharpe_train = PnL_mean_train ./ sqrt(PnL_var_train);
	finite_index_train = isfinite(PnL_sharpe_train);

	PnL_combine_train = [1:numalphas; PnL_sharpe_train']';
	PnL_combine_train = PnL_combine_train(finite_index_train, :);

	% Sort the results.
	[PnL_sort_train, order_train] = sortrows(PnL_combine_train, 2);

	% Retain the top 50 alphas.
	top_index_train = PnL_sort_train((end - num_top_alpha + 1):end, 1);
	top_sharpe_train = PnL_sort_train((end - num_top_alpha + 1):end, 2);

	% Give the top 50 alphas weight proportional to their Sharpe ratios.
	top_weight_train = top_sharpe_train / sum(top_sharpe_train);
	
	% In-sample period: backfill alpha weights to all prior days.
	% Thus strategy PnL for di <= di_start will be fully in-sample
	if di == di_start,
		alpha_weights_matrix(top_index_train, 1:di-1) = repmat(top_weight_train(:), 1, numel(1:di-1));
	end
	
	% Store weights computed on di into the di column of alpha_weights_matrix.
	alpha_weights_matrix(top_index_train, di) = top_weight_train(:);
end

% Out-of-sample period: using the weights computed on day index = di, apply them to trading on all future days.  
% Thus strategy PnL for di > di_start will be out-of-sample
for di = min(di_start + 1, data.numdates) : data.numdates,
	alpha_weights_matrix(:, di) = alpha_weights_matrix(:, di-1);
end

% TODO: using the data structure and alpha_cube structure, implement a fixed-weight strategy.
%       You may use all data and alpha_cube information up to 20131231.
%       You may pursue rolling-window strategies later, after you get experience w/ fixed-weight methods.

% ------------------------------------------------------------------------------

% The remaining steps build the preA and postA positions matrices and generate the postA's simulation result,
% which will be e-mailed to the user if the strategy finishes without errors.

% Constrain max. alpha weight.  Re-allocate weight to other alphas if needed.
% As a risk control measure, we do not allow a strategy to bet a huge amount of its capital on a single alpha.
max_weight = 0.1;
for attempt = 1:3 % or 10 times
	sum_weights = sum(abs(alpha_weights_matrix), 1);
	alpha_weights_matrix = at_nan2zero(alpha_weights_matrix ./ repmat(sum_weights, [numalphas, 1]));
	% then cap
	over_weight = (alpha_weights_matrix > max_weight);
    under_weight = (alpha_weights_matrix < -max_weight);
	alpha_weights_matrix(over_weight) = max_weight;	
    alpha_weights_matrix(under_weight) = -max_weight;	
end

% The ALL positions matrix will be generated by taking a linear weighted combination of alphas.
ALL = zeros(numstocks, numdates, 'single');

sim_start_idx = find(data.dates >= 20080101, 1, 'first');

relevantAlphasIndsCurrent = 1 : numalphas;

for di=sim_start_idx:numdates
	ALL(:,di) = (alpha_weights_matrix(:,di)' * ...
        at_nan2zero(alpha_cube.values(relevantAlphasIndsCurrent,:,di)))';
end

% Simple algorithm for in sample
all_costs = costf(ALL(:,1:di_start), data, di_start);
current = all_costs.cost;
best=current;
%the smallest weight limit.
out_limit = 10^(-7);
tol = 0.01;
    
while true

    bestLoc = 0;
    bestNeighbor = zeros(size(alpha_weights_matrix),'single');
    
    m = sum(alpha_weights_matrix(:,di_start)~=0);
    step = max_weight/m;

    for j = 1:numalphas
        if alpha_weights_matrix(j,di_start)~=0
            continue;
        end
        
        temp = alpha_weights_matrix;
        sum_weights = sum(abs(temp), 1);
        temp = at_nan2zero(temp ./ repmat(sum_weights, [numalphas, 1]))...
            .*(1-step);
        temp(abs(temp)<out_limit) = 0;         
        neighbors = repmat(temp,1,1,2);

        for k =1:2
            neighbors(j,:,k)=(-1)^k*step;
            for di=sim_start_idx:numdates
                ALL(:,di) = (neighbors(:,di,k)' * at_nan2zero(...
                    alpha_cube.values(relevantAlphasIndsCurrent,:,di)))';
            end

            all_costs = costf(ALL, data);
            cost = all_costs.cost;

            if cost<best
                best=cost;
                bestLoc = (-1)^k*j;
                bestNeighbor = neighbors(:,:,k);
            end  
        end
    end
        
    if bestLoc == 0
        break;
    end
    
    sign = bestLoc/abs(bestLoc);	
  
    cost = best;
    best = Inf;
    
    while cost < best && abs(best-cost)>tol
        best=cost;
        alpha_weights_matrix = bestNeighbor;

        temp = bestNeighbor(abs(bestLoc),:)+step*sign;
        
        if temp >max_weight
            temp =max_weight;
        elseif temp < -max_weight
            temp = -max_weight;
        end

        bestNeighbor(abs(bestLoc),:)=0;
        sum_weights = sum(abs(bestNeighbor), 1);
        bestNeighbor = at_nan2zero(bestNeighbor ./ repmat(sum_weights, ...
            [numalphas, 1])).*(1-abs(temp));
        bestNeighbor(abs(bestNeighbor)<out_limit) = 0;
        sum_weights = sum(abs(bestNeighbor), 1);
        bestNeighbor = at_nan2zero(bestNeighbor ./ repmat(sum_weights, ...
            [numalphas, 1])).*(1-abs(temp));
       
        bestNeighbor(abs(bestLoc),:)= temp;  

        for di=sim_start_idx:numdates
        	ALL(:,di) = (bestNeighbor(:,di)' * at_nan2zero(...
            	alpha_cube.values(relevantAlphasIndsCurrent,:,di)))';
        end

        all_costs = costf(ALL, data);
        cost = all_costs.cost;
    end


end          


% The preA is the positions matrix based on the combination of alphas, BEFORE imposing real-world constraints.
% The preA typically has unacceptably-high trading amounts, turnover, maximum positions, etc.
preA = ALL;

% Apply post-processing constraints.  Without the constraints below, the preA is unrealistic and not tradeable in the real world.
% (1) Position constraint #1: limit the max. abs. position of any stock to pos_limit*booksize on a given day.
% (2) Locate constraint: certain stocks are difficult to short sell.
%                        Using an estimate of total borrowable shares, assume you may borrow only locate_utilization_rate of total borrowable shares.
%                        Also, stocks whose rebate rate is below max_rebate_rate are not allowed to be traded.
% (3) Trading constraint: the strategy may only trade max_trd*vol20.
% (4) Position constraint #2: the strategy may only have a max. abs. position of max_pos*vol20.
% (5) Risk control: scale positions in factors and industries s.t. the strategy has reasonably-low net exposures each day.
% (6) Slowing: even after the steps above, turnover is generally quite high, so we smooth the desired positions over time.
%              This reduces the amount the strategy needs to trade per day, and brings turnover to a feasible level.
% Apply slippage model: trading wide-spread stocks is more costly than trading tight-spread stocks.
%                       Trades are implied from the day-over-day positions in the post-processed postA matrix, 
%                       and are then charged a fraction of the average spread per stock.
%                       Trading wide-spread illiquid names will be much more expensive than trading tight-spread liquid names.
% The resulting postA represents a realistic strategy that may be traded in the real world.
% An e-mail will be generated that plots the simulation result of the postA.

addpath('/mnt/gfs2/public/adam_metadata/helper_functions/')
[preA, postA, pre_simres, post_simres] = do_post_processing(data, preA, yesterday_preA, booksize, ...
                                                            locate_utilization_rate, max_rebate_rate, ...
															max_trd, max_pos, pos_limit, 0, ...
															delay, stockpnl, startdate, ...
															should_position_pnl, should_execution_pnl, ...
															should_slippage_pnl, should_dividend_pnl, ...
															meta_data)
															

end





% PnL simulation function
function dailypnl = calc_pnl(data, alpha,di_start)

  delay = 1;

  close = data.close;
  close(~isfinite(close)) = 0;
  alpha(~isfinite(alpha)) = 0;

  ret = close(:, (delay + 2):di_start) ./ close(:, (delay + 1):(di_start - 1)) - 1;
  ret(~isfinite(ret)) = 0;
  ret_seq = sum(alpha(:, 1:(di_start - (delay + 1))) .* ret, 1);
  dailypnl = [transpose(zeros((delay + 1), 1)), ret_seq];
  dailypnl = dailypnl(:);

end

function cumpnl = calcCumPnl(dailypnl)
    %dailypnl(numAlphas,days);
    
    cumpnl = dailypnl;
    for n = 2:size(cumpnl,2)
        cumpnl(:,n) = cumpnl(:,n) + cumpnl(:,n-1);
    end
    
end

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

end


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

function alpha = calcScores(dailypnl)
% scoring system
% factors: fs
% weights: ws

    alpha.sharpeR =  calcSharpeR(dailypnl );
    alpha.cumpnl = calcCumPnl(dailypnl);
	alpha.maxDD =  calcMaxDD( alpha.cumpnl );
	alpha.kRatio = calcKratio(alpha.cumpnl,252); 
    
    ws = [10,-0.1,0.1];
	alpha.score = ws*[alpha.sharpeR;alpha.maxDD;alpha.kRatio];

end


function sharpeR = calcSharpeR(dailypnl)

    dailypnl(~isfinite(dailypnl)) = 0;
    %dailypnl(numAlphas,days);
    
    % dolar neutral so no substraction of risk-free return.
    sharpeR = mean(dailypnl,2)./std(dailypnl,0,2);
    sharpeR(~isfinite(sharpeR)) = 0;
end

function alpha = costf(s_postA, data, di)

%     tradable cost
      C = 1*10^9;
      %upper-bound
      ub = sum(abs(s_postA),1)- C;
      ub(ub<0) = 0;
      alpha.ubCost = sum(ub(:));
      %clear ub;
      
      delay = 1;          
      %Limited trading
      shares=(s_postA(:,2:end-delay)-s_postA(:,1:end-1-delay))./data.close(:,2+delay:di);
      shares(~isfinite(shares))=0;
      
      %round trades
      shares = round(shares./100).*100;      

      %volume-based trading limit
      t = 0.01;      
      lt = abs(shares(:,21-1-delay:end)) -  t.*data.vol20(:,21:di);
      lt(lt<0)=0;
      alpha.ltCost = sum(lt(:));   
      %clear lt;

      %limited short sale
      %utilization rate
      u = 0.1;
      ss =  shares;
      ss(s_postA(:,1:end-1-delay)>0)=0;
      ss(ss>0)=0;
      lss = abs(ss) - u*data.locates_ms_avshs(:,2+delay:di);
      lss(lss<0) = 0;
      alpha.lssCost = sum(lss(:));  
      clear lss;
      %clear ss;

      %No hard to borrow trading
      nhb =  shares;
      nhb(data.locates_ms_rebate(:,2+delay:di)>=-10)=0;
      alpha.nhbCost = sum(abs(nhb(:)));
      %clear nhb;
      
      %position limit
      %p is captial based position limit
      p = 0.5/100;
      pl1 = abs(s_postA)-p*C;
      pl1(pl1<0)=0;
      alpha.pl1Cost = sum(pl1(:));
      %clear pl1;
      
      %m volume based position limit
      m = 0.15;
      pl2 = abs(s_postA(:,21-delay:end-delay))-m.*data.vol20(:,21:di).*data.close(:,21:di);
      pl2(pl2<0) = 0;
      alpha.pl2Cost = sum(pl2(:));
      %clear pl2;
      
      %Dollar neutrality
      epsilon_dollar = 0.05*C;
      dn = abs(sum(s_postA)) - epsilon_dollar;
      dn(dn<0)=0;
      alpha.dnCost = sum(dn(:));
      %clear dn;
      
      %industry neutrality
      epsilon_industry = 0.05*C;   
    
%       indu = unique(data.industry(:));
%       induCost = 0;
%     
%         for i = 2:size(indu,1)
%             i_postA = s_postA;
%             i_postA(data.industry~=indu(i))=0;
%             in = abs(sum(i_postA))-epsilon_industry;
%             in(in<0) = 0;
%             induCost = induCost + sum(in,2);
%         end
    induCost = 0;    
    for i = 1:size(data.industry,2)
        i_postA = s_postA;
        i_postA(~data.industry(i).max)=0;
        in = abs(sum(i_postA))-epsilon_industry;
        in(in<0) = 0;
        induCost = induCost + sum(in,2);
    end
    alpha.induCost = induCost;
 
    %Q: 0 for industry?  
    
    %factor neutrality? TODO    

    alpha.cost_tradable = alpha.ubCost + alpha.ltCost +  alpha.lssCost...
      + alpha.nhbCost + alpha.pl1Cost ...
      + alpha.pl2Cost + alpha.dnCost + alpha.induCost;
  
  
    dailypnl = calc_pnl(data, s_postA, di_start);   
       
    scores = calcScores(dailypnl');
    alpha.score = scores.score;
    alpha.sharpeR =  scores.sharpeR;
	alpha.maxDD =  scores.maxDD;
	alpha.kRatio = scores.kRatio;
    
% don't include performance cost for now.
    %score =0;
    alpha.cost_performance = -alpha.score;    
    alpha.cost = alpha.cost_tradable + alpha.cost_performance;
end

