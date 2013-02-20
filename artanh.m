function out=artanh(x)
%artanh function
%abs(x)<1

if abs(x)>=1
    fprintf('abs of input param must be less than 1\n');
    out=[];
else
    out=0.5*log((1+x)/(1-x));
end

