function [xs,ys] = shrink2(x,y,lambda)
%--------------------------------------------------------------------------
%   This program is computing the generalized shrinkage (2D).
%    
%   Usage: [xs,ys] = shrink2(x,y,lambda);
%
%   Inputs: 
%       - x: 2d data
%       - y: 2d data
%       - lambda: threshold used to do shrinkage
%
%   Outputs: 
%       - xs: shrinkage of x
%       - ys: shrinkage of y
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------



[xLen,yLen]=size(lambda);
if xLen==1
    % isotropic 
    s = sqrt(x.*conj(x)+y.*conj(y));
    ss = s-lambda;
    ss = ss.*(ss>0);

    s = s+(s<lambda);
    ss = ss./s;

    xs = ss.*x;
    ys = ss.*y;
else
    i=sqrt(-1);
    s=x+i*y;
    ss=wthresh(s,'s',lambda);
    xs=real(ss);
    ys=imag(ss);
end