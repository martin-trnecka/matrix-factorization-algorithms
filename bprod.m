function I = bprod(A, B)
% BPROD a Boolean product of two matrices
% usage: A = bprod(B, C) returns I = A \circ B.

I = min(1, double(A)*double(B));

end
