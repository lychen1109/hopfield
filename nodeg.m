function out=nodeg(u,u0)
%node output

out=0.5*(1+tanh(u/u0));
