function [ A, B ] = GreConDPlus( M, w, no_of_factors )
% GRECONDPLUS implements GreConD+ algorithm for Boolean matrix factorization 

% usage: [A, B] = GreConD(I, w, no_of_factors);
% returns A \approx B = I
% w drives a penalization of overcover error

% if you are using this implementation please cite the following work
% Belohlavek R., Trnecka M.: 
% A new algorithm for Boolean matrix factorization which admits overcovering
% Discrete Applied Mathematics 249(20)(2018), 36-52.

M = logical(M); % logical type is required
[m, n] = size(M);
coverage = zeros(m, n);
U = M;
factors = 0;

aa = logical([]);
bb = logical([]);
E = logical([]); % row expansion of factors
F = logical([]); % col expansion of factors

% spusteni vypoctu
while any(any(U))
 v = 0;
    d = false(1,n);
    d_old = false(1,n);
    d_mid = false(1,n);
    e = true(m,1); % extent for speed closure
    
    while 1
        for j=1:n
            if ~d(j)
                % computes the value of the cover function for the candidate factor
                % inline function for speed
                % arrow down (speed version)
                a = e & M(:,j);
                % arrow up
                sum_a = sum(a);
                if sum_a*n > v % check the size of upper bound
                    b = all(M(a,:),1);
                
                    if sum_a*sum(b) > v % check the size of upper bound
                        cost = sum(sum(U(a,b)));
                
                        if cost > v
                            v = cost;
                            d_mid = b;
                            c = a;
                        end
                    end
                end
            end
        end
        
        d = d_mid;
        e = c;
        
        if all(d==d_old)
            break;
        else
            d_old = d;
        end
    end
    
    % [e, f] is expansion of [C, D]
    [e, f] = expansion(c, d, U, M, w);
    c = or(c, e);
    d = or(d, f);
    
    aa = [aa, c];
    bb = [bb; d];
    E = [E, e];
    F = [F; f];
    
    factors = factors + 1;
    display(factors);
    
    % end if the no. of factors is reached
    if nargin==3 && factors==no_of_factors
        break;
    end
    
    % odebrani z U a pridani do cover matice
    U(c,d) = 0;
    coverage(c,d) = coverage(c,d) + 1;
    
    
    % check if there is factor than can be removed
    for k=1:size(bb, 1)
        if k >= size(bb, 1)
            break;
        end
        
        c = coverage(aa(:,k), bb(k,:));
        
        if min(min(c(M(aa(:,k), bb(k,:))))) >= 2
            disp('factor was removed');
            
            aa(:,k) = [];
            bb(k,:) = [];
            E(:,k) = []; % del expansion
            F(k,:) = [];
            coverage(aa(:,k), bb(k,:)) = coverage(aa(:,k), bb(k,:)) - 1;
        end
    end
    
    % check if the overcover error may decrease
    for k=1:size(bb, 1)
        % loop over all col in expansion
        for j=1:n
            if F(k,j)
                z = false(1,n);
                z(j) = 1;
                c = coverage(aa(:,k), z);

                if min(min(c(M(aa(:,k), z)))) >= 2
                    disp('col expansion was removed');
                    b(k,j) = 0;
                    coverage(aa(:,k), z) = coverage(aa(:,k), z) - 1;
                    F(k,j) = 0;
                end
            end
        end
        
        % loop over all row in expansion
        for j=1:m
            if E(j,k)
                z = false(m,1);
                z(j) = 1;
                c = coverage(z, bb(k,:));

                if min(min(c(M(z, bb(k,:))))) >= 2
                    disp('row expansion was removed');
                    aa(j,k) = 0;
                    coverage(z, bb(k,:)) = coverage(z, bb(k,:)) - 1;
                    E(j,k) = 0;
                end
            end
        end
    end
      
end


A = logical(aa);
B = logical(bb);

end



% expansion of factor
function [E, F] = expansion(a, b, U, M, w)

[m, n] = size(U);
E = false(m,1);
F = false(1,n);

while 1
    cover = sum(sum(U(a, b)));
    overcover = sum(sum(~M(a, b)));
    cost = cover - w * overcover;
    
    co = 0; % control if the row or colum will be expanded
    
    % cost of columns that are not in B
    price = cover + sum(U(a,:)) - w * (overcover + sum(~M(a,:),1));
    price(b) = -Inf;
    [~,col] = max(price);
   
    %  cost of the best col
    if cost < price(col)
        co = 2;
    end
        
    % cost of rows that are not in A
    price = cover + sum(U(:,b),2) - w * (overcover + sum(~M(:,b),2));
    price(a) = -Inf;
    [~,row] = max(price);
   
    % cost of the best row
    if cost < price(row)
        co = 1;
    end
    
    if co == 0 % no expansion
        break;
    elseif co == 1 % add row
        a(row) = 1;
        E(row) = 1;
    elseif co == 2 % add col
        b(col) = 1;
        F(col) = 1;
    end
end

end