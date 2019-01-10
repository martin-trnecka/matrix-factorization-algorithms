% example of usage
load mushroom.mat

tic();
[A,B] = GreConD(mushroom);
toc();

I = bprod(A,B); 