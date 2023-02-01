function [im,u] = matlab_test()

A = [1,2,3;4,5,6;7,8,9];
Fx(A)
Fy(A)
Bx(A)
By(A)

img_name = 'ultrasound.jpg';
clean = imread(img_name);
clean = im2double(clean);

    sig = 15; %noise level
    %add noise:
    im = imnoise(clean,'gaussian',0,(sig/255)^2);

    lambda1 = 5;
    mu1 = 0;
    sigma1 = 2;
    xi1 = 0*1000;
    lambda2 = 20;
    mu2 = 1;
    sigma2 = 2;
    xi2 = 1000;
    lambda3 = 20;
    mu3 = 1;
    sigma3 = 2;
    xi3 = 1000;

    %game start:
    [n,m] = size(im);
    u = im;
    wx = zeros(n,m); wy = wx; bx = wx; by = wx;
    v = wx;
    w2x = wx;
    w2y = wx;
    b2x = wx;
    b2y = wx;
    z = wx;
    w3x = wx;
    w3y = wx;
    b3x = wx;
    b3y = wx;

    maxit = 300;

    lhs_u = lambda1 + 4*sigma1;
E=[];
    for k =1:maxit

        for k2 = 1:5
        rhs = lambda1*im - sigma1*Bx(wx-bx)- sigma1*By(wy-by)+ sigma1*(Bx(Fx(u)) + By(Fy(u)) + 4*u);
        u = rhs./lhs_u;
        end

        temp1 = Fx(u) + bx; temp2 = Fy(u) + by;
        [wx,wy] = shrink2b(temp1,temp2,1/sigma1);

        bx = bx + Fx(u) - wx;
        by = by + Fy(u) - wy;
    
        E1 = sum(sum((im-u).^2));
        TV = Fx(u).^2 + Fy(u).^2;
        TV = sum(sum(sqrt(TV)));
        E = [E; E1+TV];
    end
    figure; plot(E);



end

function [wx,wy] = shrink2b(sx,sy,sig)
s = sqrt(sx.^2 + sy.^2 + 1e-9);
t = max(s-sig,0);
wx = t.*(sx./s);
wy = t.*(sy./s);


end

% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1])-u;
end

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0])-u;
end

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u - circshift(u,[0 1]);
end

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u - circshift(u,[1 0]);
end