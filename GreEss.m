function [ A, B ] = GreEss( I )
% GreEss algorithm
% The algorithm for a given binary matrix I finds the Boolean decomposition of I into matrices A and B.
% namely the min(min(min(1, double(A)*double(B))==I)) is equal to 1

% If you use this code, please cite
% Belohlavek R., Trnecka M.:
% From-Below Approximations in Boolean Matrix Factorization: Geometry and New Algorithm.
% Journal of Computer and System Sciences 81(8)(2015),
% 1678�1697, DOI 10.1016/j.jcss.2015.06.002.
% note, the following implementation is significantly optimized against the implementation used in the paper
% the original implementation can be found on http://trnecka.inf.upol.cz/publications/implementations/GreEss.zip

% usage:
% [A, B] = GreEss(I)
% Matrix I should be a MATLAB logical type

% auxiliary methods in C:
% c_ess.c  ... computes the Ess matrix
% source code in C can be compiled via mex function in MATLAB

% convert to the logical type
I = logical(I);
U = I; % U represents the uncovered part of I
[m, n] = size(I);
empty = false(m,n);

% Ess represent the essential part of I
disp('Compute Ess part...');
Ess = c_ess(I);
disp('Done');

% compute intervals (each interval includes candidates for factor)
disp('Compute candidate intervals...');
[C, D] = computeIntervals(Ess, I);

no_of_intervals = size(C, 2);
exclusion = false(1, no_of_intervals);
number_of_factors = 0;
A = logical([]);
B = logical([]);

% run until all elements of I are covered by some factor
while any(any(U))

    area = 0;

    for k=1:no_of_intervals
        % take only one factor from each ess interval
        if exclusion(k)
            continue;
        end

        % compute interval [c, d]
        d = all(I(C(:,k),:), 1);
        c = all(I(:,D(k,:)), 2);

        % create area which represent interval [c, d] in I
        submatrix = empty;
        submatrix(c, d) = 1;
        submatrix = submatrix & I;

        % compute best factor (simpified via GreConD algorithm)
        [E, F] = get_factor(submatrix, U);

        %if (min(x==E) == 0) || (min(y==F) == 0)
        %    break;
        %end

        % compute the number of nonzero covered entries
        p = sum(sum(U(E,F)));

        if p > area
            area = p;
            vE = E;
            vF = F;
            selectedInterval = k;
        end
    end

    number_of_factors = number_of_factors + 1;
    display(number_of_factors);

    A = [A, vE];
    B = [B; vF];

    % remove phase
    U(vE, vF) = 0;
    exclusion(selectedInterval) = 1;
end
end



% GreConD algorithm for navigation in particular interval
function [c, d] = get_factor(submatrix, U)

[m, n] = size(submatrix);
d_old = false(1,n);
d_temp = false(1,n);
d = false(1,n);
c = false(m,1);
e = true(m,1); % extent for speed closure
v = 0;

atr = find(sum(submatrix)>0); % only not covered attributes

%run until there are changes
while 1
    for j=atr
        if ~d(j)
            d(j) = 1;

            a = e & submatrix(:,j);
            sum_a = sum(a);

            if sum_a*n > v
                b = all(submatrix(a,:), 1);
                if sum_a*sum(b) > v
                    area = sum(sum(U(a,b)));
                    d(j) = 0;

                    if(area > v)
                        v = area;
                        d_temp = b;
                        c = a;
                    end
                end
            end

            d(j) = 0;
        end
    end

    d = d_temp;

    if(all(d==d_old))
        break
    else
        d_old = d;
    end
end

end








% computes all intervals
function [ A, B ] = computeIntervals( Ess, I )

[~, n] = size(I);
U = Ess;
A = logical([]);
B = logical([]);

while any(any(U))
    v = 0;
    d = false(1,n);
    d_old = false(1,n);
    d_mid = false(1,n);
    finalRemoveA = [];
    finalRemoveB = [];

    atr = find(sum(U)>0); % only not covered attributes

    while 1
        for j=atr
            if ~d(j)
                d(j)=1;

                a = all(Ess(:,d), 2);
                sum_a = sum(a);

                if sum_a*n > v % simple speedup
                    b = all(Ess(a,:), 1);
                    if sum_a*sum(b)>v % simple speedup
                        removeB = all(I(a,:),1);
                        removeA = all(I(:,removeB),2);

                        area = sum(sum(U(a,b)));

                        d(j)=0;

                        if(area > v)
                            v = area;
                            d_mid = b;
                            c = a;
                            finalRemoveA = removeA;
                            finalRemoveB = removeB;
                        end
                    end
                end

                d(j)=0;
            end
        end

        d = d_mid;

        if all(d==d_old)
            break;
        else
            d_old = d;
        end

    end

    A = [A, c];
    B = [B; d];

    % remove covered elements
    U(finalRemoveA, finalRemoveB) = 0;
end

end
