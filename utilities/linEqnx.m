function [x, y] = linEqnx(f1, f2)

%Finds intersection of functions of type
%y=mx+b or x=c

m1 = f1(1);
b1 = f1(2);
c1 = f1(3);
m2 = f2(1);
b2 = f2(2);
c2 = f2(3);


if isnan(c1) && isnan(c2)
    
    if m1==m2
        x = NaN;
        y = NaN;
    else
    
    x = (b2-b1)/(m1-m2);
    y = m1*x + b1;
    end
    

elseif ~isnan(c1) && isnan(c2)
    x = c1;
    y = (m2*c1)+b2;
    
elseif isnan(c1) && ~isnan(c2)
    x = c2;
    y = (m1*c2)+b1;
    
else
    x = NaN;
    y = NaN;
end


end
    
    





