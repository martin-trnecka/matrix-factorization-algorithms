% example of usage
load mushroom.mat

% GreConD algorithm
tic();
[A,B] = GreConD(mushroom);
toc();

I = bprod(A,B);


% GreEss algorithm
tic();
[C,D] = GreConD(mushroom);
toc();

I = bprod(C,D);
