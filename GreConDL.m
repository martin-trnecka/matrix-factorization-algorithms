function [ A, B ] = GreConDL( I, L )
% GRECONDL implements algorithm GreConD for L-set, use
% Lukasiewicz conjunction and residum
% I ... matrix with ordinal data (from scale L)
% L ... scale, e.g. [0, 0.25, 0.5, 0.75, 1]

% if you are using this implementation please cite the following work
% Belohlavek R., Vychodil V.: 
% Factor analysis of incidence data via novel decomposition of matrices. 
% Lecture Notes in Artificial Intelligence 5548(2009), 83–97. 
% [Springer, ISSN 0302-9743, DOI 10.1007/978–3–642–01815–2_8]

[~, n] = size(I);
number_of_degrees = size(L, 2);
U = (I>0); % still uncovered elements
A = [];
B = [];
factors = 0;

% loop until all nozero elements are covered
while any(any(U))
    
    D2 = zeros(1, n);
    V = -Inf;
    v_max = 0;
    j_max = 0;
    a_max = 0;
    
    atr = find(sum(U)>0); % only not covered attributes
    
    % loop until the cover is changing
    while v_max > V
        V = v_max;
        C = arrowDown(D2, I);
        D = arrowUp(C, I);
        
        % loop over all columns
        for j=atr
            % loop over all degrees in L
            for a=1:number_of_degrees
                
                D1 = D;
                if(L(a) > D1(j))
                    D1(j) = L(a);
                    C = arrowDown(D1, I);
                    D1 = arrowUp(C, I);
                    
                    % compute the cost function
                    luk = max(C + D1 - 1, 0);
                    velikost = sum(sum(U(luk >= I)));
                    
                    % maximize the cost function
                    if(velikost > v_max)
                        v_max=velikost;
                        j_max=j;
                        a_max=a;
                    end
                end
            end
        end
        
        D2 = D; % store actual intent of factor
        D2(j_max) = L(a_max);
    end
    
    % store computed factor
    C = arrowDown(D, I);
    A = [A, C];
    B = [B; D];
    
    factors = factors + 1;
    disp(factors);
    
    % remove elements covered by [C, D]
    luk = max(C + D - 1, 0);
    rozdil = abs(I - luk);
    U(rozdil <= 0.000001) = 0;
    
end
end


% arrow down (Lukasiewicz)
function [ D ] = arrowDown( C, matice )
%[m, ~] = size(matice);
%factor = min(1-repmat(C, m, 1) + matice, 1);

%     for i=1:m
%         for j=1:n
%             faktor(i,j)=min(1-C(j)+matice(i,j),1);
%         end
%     end
factor = min(1-C+matice,1);
D = min(factor, [], 2);
end


% arrow up (Lukasiewicz)
function [ D ] = arrowUp( C, matice )
%[~, n]=size(matice);
%factor = min(1-repmat(C, 1, n) + matice, 1);

%     for i=1:m
%         for j=1:n
%             faktor(i,j)=min(1-C(i)+matice(i,j),1);
%         end
%     end

factor = min(1-C+matice,1);
D = min(factor, [], 1);
end



