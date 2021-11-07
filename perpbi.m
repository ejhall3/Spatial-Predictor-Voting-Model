function [y, mperp, b, c] = perpbi(p1, p2)

[~, col] = size(p1);

for i =1:col
    mid(i) = (p1(i)+p2(i))/2;
end

p1x = p1(1);
p1y = p1(2);
p2x = p2(1);
p2y = p2(2);

m = (p2y - p1y) / (p2x - p1x);

if m == 0
    c = mean([p1x, p2x]);
else
    c = NaN;
end
    

mperp = -1/m;

x = linspace(-5, 5, 1000);
b = (-1*mperp*mid(1)) + mid(2);
y = (mperp*x) + b;

end
