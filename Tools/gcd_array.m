function v = gcd_array(v)
% finds the greatest common divisor of all elements in a vector v

arguments
    v (1,:) double
end

while length(v) > 1
    v = [gcd(v(1), v(2)) v(3:end)];
end

end