function y = two_cmpl (x)
y = bin_add(1-x,[zeros(1,numel(x)-1) 1]);
end