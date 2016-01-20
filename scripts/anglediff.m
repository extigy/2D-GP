function d = anglediff(th1, th2)
        if nargin < 2
            d = th1;
        else
            d = th1 - th2;
        end
        d = mod(d+pi, 2*pi) - pi;
    end