function [ret] = MOD(x,y)
    if x < 0
        ret = x+y;
    else
        ret = mod(x,y);
    end 
end
