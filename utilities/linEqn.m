function [out] = linEqn(m, b, c)

%Return y(x)= mx+b evaluated at x=c

if isnan(c)
    out = NaN;
else  
    out = (m*c)+b;
end

end


