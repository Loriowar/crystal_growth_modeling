function [ ret ] = is_integer( a )
if a == round(a)
    ret = true;
else
    ret = false;
end
end

