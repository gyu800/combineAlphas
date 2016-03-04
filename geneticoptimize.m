function sol = geneticoptimize( alphas, data, costf, maxiter, popsize, ...
   mutprob, elite, step, tol)
    
    numAlpha = size(alphas,1);
    %build initial pop
    for i = 1:popsize
        pop(i,1).sol = rand(numAlpha, 1)>0.5;
        pop(i,1).ws = zeros(numAlpha,1);
        pop(i,1).ws(pop(i,1).sol)=1;
        pop(i,1).ws(pop(i,1).sol)=1/sum(pop(i,1).ws);
        pop(i,1).cost = [];
        pop(i,1).alpha = [];
    end
    
    
    topelite = elite*popsize;
    
    for i = 1:maxiter
        display(datestr(now));
        fprintf('genetic op iter %d\n',i);
        
        for j = 1:popsize
            if isempty(pop(j).cost)
                display(datestr(now));
                fprintf('calc %d pop cost\n',j);
                
                results = hillclimb(alphas(pop(j).sol),pop(j).ws(pop(j).sol)...
                    , data,costf, step, maxiter, tol);
                pop(j).cost = results.alpha.cost;
                pop(j).alpha = results.alpha;
                pop(j).ws(pop(j).sol) = results.ws;
                pop(j).sol(pop(j).ws==0) = false;
            end
        end
        [sortedPop,sortIndex] = sort([pop(:).cost]);
        pop = pop(sortIndex);
        save(['popp',num2str(i),'.mat'],'pop');
        display(datestr(now));
        fprintf('pop1.cost=%d\n',pop(1).cost);
        for j = topelite+1:popsize
            if rand<mutprob
                %mutation
                mutated = pop(randi(topelite));
                k = randi(numAlpha);
                if mutated.sol(k)
                    mutated.ws(k) = 0;
                else
                    mutated.ws(k) = step*2;
                end
                mutated.sol(k)=~mutated.sol(k); 
                mutated.ws = mutated.ws./sum(mutated.ws);
                mutated.cost = [];
                pop(j) = mutated;
            else
                %Crossover
                c1 = pop(randi(topelite));
                c2 = pop(randi(topelite));
                k = randi(numAlpha-1);
                c1.sol(k+1:end)=c2.sol(k+1:end);
                c1.ws(k+1:end)=c2.ws(k+1:end);
                c1.ws = c1.ws./sum(c1.ws);
                c1.cost = [];
                pop(j) = c1;
            end
        end
                
        
    end
    
     sol = pop;

end