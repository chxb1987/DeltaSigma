function y = sgn_xtnd (x)
y = [x(1) x(1:numel(x)-1)];
end