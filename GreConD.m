function [ A, B ] = GreConD( I, no_of_factors )
% GRECOND implements GreConD algorithm for Boolean matrix factorization 

% usage: [A, B] = GreConD(I);
% returns A \circ B = I (if the no. of factors is not limited)

% if you are using this implementation please cite the following work
% Belohlavek R., Vychodil V.: 
% Discovery of optimal factors in binary data via a novel method of matrix decomposition. 
% Journal of Computer and System Sciences 76(1)(2010), 3-20. 

M = logical(I);
[m, n] = size(M);
U = M;
k = 0;

A = logical([]);
B = logical([]);

while any(any(U))
    v = 0;
    d = false(1,n);
    d_old = false(1,n);
    d_mid = false(1,n);
    e = true(m,1); % extent for speed closure
    
    atr = find(sum(U)>0); % only not covered attributes
    
    while 1
        for j=atr
            if ~d(j)
                % computes the value of the cover function for the candidate factor
                % inline function for speed
                % arrow down (speed version)
                a = e & M(:,j);
                % arrow ups
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
    
    A = [A, c];
    B = [B; d];
    
    k = k + 1;
    display(k);
    
    % end if the no. of factors is reached
    if nargin==2 && k==no_of_factors
        break;
    end
    
    % delete already covered part
    U(c, d) = 0;
end
end