function z = bin_add (x, y)
if (numel(x) ~= numel(y))
    if(numel(x) < numel(y))
        x = [zeros(1,numel(y)-numel(x)) x];
    else
        y = [zeros(1,numel(x)-numel(y)) y];
    end
end
z = zeros(1,numel(x));
for i=numel(x):-1:2
    z(i-1) = mod(z(i-1)+floor((z(i)+x(i)+y(i))/2),2);
    z(i) = mod(z(i)+x(i)+y(i),2);
end
z(1) = mod(z(1)+x(1)+y(1),2);
end