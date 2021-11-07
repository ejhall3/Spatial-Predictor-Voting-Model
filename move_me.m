function v = move_me(v,a)
if nargin <2 
    a = 0;
end
 v = [v(v ~= a), v(v == a)];
end