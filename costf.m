function alpha = costf(s_postA, data)
%     display(datestr(now));
%     display('cost calculation');
% tradable cost
      C = 1*10^9;
      %upper-bound
      ub = sum(abs(s_postA),1)- C;
      ub(ub<0) = 0;
      alpha.ubCost = sum(ub(:));
      %clear ub;
      
      delay = 1;          
      %Limited trading
      shares=(s_postA(:,2:end-delay)-s_postA(:,1:end-1-delay))./data.close(:,2+delay:end);
      shares(~isfinite(shares))=0;
      
      %round trades
      shares = round(shares./100).*100;      

      %volume-based trading limit
      t = 0.01;      
      lt = abs(shares(:,21-1-delay:end)) -  t.*data.vol20(:,21:end);
      lt(lt<0)=0;
      alpha.ltCost = sum(lt(:));   
      %clear lt;

      %limited short sale
      %utilization rate
      u = 0.1;
      ss =  shares;
      ss(s_postA(:,1:end-1-delay)>0)=0;
      ss(ss>0)=0;
      lss = abs(ss) - u*data.locates_ms_avshs(:,2+delay:end);
      lss(lss<0) = 0;
      alpha.lssCost = sum(lss(:));  
      clear lss;
      %clear ss;

      %No hard to borrow trading
      nhb =  shares;
      nhb(data.locates_ms_rebate(:,2+delay:end)>=-10)=0;
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
      pl2 = abs(s_postA(:,21-delay:end-delay))-m.*data.vol20(:,21:end).*data.close(:,21:end);
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
    for i = 1:size(data.indu,2)
        i_postA = s_postA;
        i_postA(~data.indu(i).max)=0;
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
  
  
    dailypnl = calc_pnl(data, s_postA);   
       
    scores = calcScores(dailypnl', data.ws);
    alpha.score = scores.score;
    alpha.sharpeR =  scores.sharpeR;
	alpha.maxDD =  scores.maxDD;
	alpha.kRatio = scores.kRatio;
    
% don't include performance cost for now.
    %score =0;
    alpha.cost_performance = -alpha.score;    
    alpha.cost = alpha.cost_tradable + alpha.cost_performance;
end

