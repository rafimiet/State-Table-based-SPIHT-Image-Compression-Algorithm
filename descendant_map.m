function [m,LB] = descendant_map(x,y,S,T)
    % x,y are the indices of current block
    % S is the input image
    % m = 1 if descendants significant
    % LB is output contains block addresses of descendants which are BL < 7

[r,c] = size(S);
xx = x;yy = y;i = 1;m = 0;
while ((xx <= r/4) && (yy(1) <= c/4))
    mx = max(max(abs(S(4*(2^(i-1))*x-4*(2^(i-1))+1:4*(2^(i-1))*x, 4*(2^(i-1))*y-4*(2^(i-1))+1:4*(2^(i-1))*y))));
    if mx >= T
        m = 1;
        break
    end
    xx = 2*xx; yy = 2*yy;
    i = i + 1;
end
LB = [2*x - 1, 2*y - 1;2*x - 1, 2*y;2*x, 2*y - 1;2*x, 2*y];