function solution = simplex( alphas, ws, data, costf, maxW, maxiter, tol )
      
    n = size(alphas,1);    
    alphaNames = [alphas(:).name]';
    
    % See initial result is
    s_postA = combineAlphas(alphaNames(ws~=0),ws(ws~=0), data);
    strategy = costf(s_postA, data);
    current = strategy.cost;

    best=current;
    best_strategy = strategy;    
    i = 0;    
    while maxiter>0
        %maxiter = maxiter - 1;
        display(datestr(now));
        i = i+1;
        fprintf('simplex iter = %d\n',i);
        
        bestLoc = 0;
        m = sum(ws~=0);
        step = maxW/m;
        
        for j = 1:n
%             if ws(j)~=0
%                 continue;
%             end
%             display(datestr(now));
%             fprintf('calc %d neighbor cost\n',j);  
            temp = ws;
            temp=temp./sum(abs(temp)).*(1-step);
            temp(~isfinite(temp)) = 0;
            %temp(abs(temp)<10^(-6)) = 0;
            neighbors = repmat(temp,1,2); 
 
            for k =1:2
                neighbors(j,k)=(-1)^k*step;
                s_postA = combineAlphas(alphaNames(neighbors(:,k)~=0),neighbors(neighbors(:,k)~=0,k),data);
                strategy = costf(s_postA, data);
                cost = strategy.cost;

                if cost<best
                    best=cost;
                    best_strategy = strategy;
                    bestLoc = (-1)^k*j;
                    fprintf('simplex %d negibors cost=%d\n',bestLoc, best);
                    %save('loc.mat','bestLoc');
                    save('best2.mat','best');
                    save('neighbors2.mat','neighbors');
                end  
            end
        end
        
        if bestLoc == 0
            break;
        end
        sign = bestLoc/abs(bestLoc);
        ws = ws./sum(abs(ws)).*(1-step);
        ws(~isfinite(ws)) = 0;
        ws(abs(bestLoc))=step*sign;

        neighbor = ws;     
        cost = best;
        best = Inf;
        strategy = best_strategy;
        
        while cost < best && abs(best-cost)>tol
            best=cost;
            best_strategy = strategy;
            ws = neighbor;
            save('best2.mat','best');
            save('ws2.mat','ws');
            
            temp = neighbor(abs(bestLoc))+step*sign;
            if temp >maxW
                temp =maxW;
            elseif temp < -maxW
                temp = -maxW;
            end
            neighbor(abs(bestLoc))=0;
            neighbor=neighbor./sum(abs(neighbor)).*(1-abs(temp)); 
            neighbor(~isfinite(neighbor)) = 0;
            neighbor(abs(neighbor)<10^(-6)) = 0;
            neighbor=neighbor./sum(abs(neighbor)).*(1-abs(temp)); 
            neighbor(~isfinite(neighbor)) = 0;
            neighbor(abs(bestLoc))= temp;  
            %neighbor(~isfinite(neighbor)) = 0;

            s_postA = combineAlphas(alphaNames(neighbor~=0),neighbor(neighbor~=0),data);
            strategy = costf(s_postA, data);
            cost = strategy.cost;
        end
        
%         m = sum(ws~=0);
%         step = maxW/m;
%         for j = 1:n
%             if ws(j)==0
%                 continue;
%             end
%             
%             temp = ws;
%             temp=temp./sum(abs(temp)).*(1-step);
%             temp(~isfinite(temp)) = 0;
%             %temp(abs(temp)<10^(-6)) = 0;
%             neighbors = repmat(temp,1,2); 
%  
%             for k =1:2
%                 neighbors(j,k)=neighbors(j,k)+(-1)^k*step;
%                 s_postA = combineAlphas(alphaNames(neighbors(:,k)~=0),neighbors(neighbors(:,k)~=0,k),data);
%                 strategy = costf(s_postA, data);
%                 cost = strategy.cost;
% 
%                 if cost<best
%                     best=cost;
%                     best_strategy = strategy;
%                     bestLoc = (-1)^k*j;
%                     fprintf('simplex %d negibors cost=%d\n',bestLoc, best);
%                     save('loc.mat','bestLoc');
%                     save('best.mat','best');
%                     save('neighbors.mat','neighbors');
%                 end  
%             end
%             
%         end



    end           
  
    
    solution.ws = ws;
    solution.alpha = best_strategy;
    
end

