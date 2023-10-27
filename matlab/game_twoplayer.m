function [u,v,z]=game_twoplayer(Img,lambda1,mu1,lambda2,mu2,sigma1,sigma2,xi1,xi2)
%--------------------------------------------------------------------------
% v is Gx
% z is Gy


%% switches
update_dynamically=1;
%%
% initialization
[n,m] = size(Img);
u=Img;wx=0*Img;wy=wx;bx=wx;by=wx;
v=zeros(n,m); w2=v; b2=v;
z=zeros(n,m); 
% Forward difference matrices
[xLen,yLen]=size(Img);
[DDx,DDy] = first_order_differences(xLen,yLen);

% Build Kernels: use the fft algorithm: (5 point stencil)
uker = 0*Img;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;  

% compute the fft of the left hand side 

uker1 = lambda1+(mu1+sigma1 + xi1)*fft2(uker);
uker2 = lambda2 + xi2 +(mu2+sigma2)*fft2(uker);


Imgx = DDx*Img; Imgy = Img*DDy';
[Imgy,Imgx] = gradient(Img);

%%
for i=1:300

    u0=u;
    v0=v; z0=z;
    rhs1 = lambda1.*Img+sigma1*DDx'*(wx-bx)+sigma1*(wy-by)*DDy + xi1*DDx'*v + xi1*z*DDy;
    newu = ifft2(fft2(rhs1)./uker1);
    if update_dynamically == 1
    u = newu;
    end
    
    rhs2 = lambda2.*Imgx + xi2.*DDx*u +sigma2*DDx'*(w2-b2);
    newv = ifft2(fft2(rhs2)./uker2);
    if update_dynamically == 1
    v = newv;
    end
    
    
    rhs3 = lambda2.*Imgy + xi2*u*DDy' +sigma2*(w2-b2)*DDy;
    newz = ifft2(fft2(rhs3)./uker2);
    if update_dynamically == 1
        z = newz;
    end    
    

    
    
    
    err1=norm(newu-u0,'fro')/norm(u,'fro');
    err2=norm(newv-v0,'fro')/norm(u,'fro');
    err3=norm(newz-z0,'fro')/norm(u,'fro');
    err = err1+err2+err3;
    if mod(i,10)==0
        disp(['iterations: ' num2str(i) '!  ' 'error is:   ' num2str(err1+err2+err3)]);
    end
    
    % check the stopping criterion
    if err<10^(-4)
        break;
    end
    
    % w problem
    temp1=DDx*u+bx;temp2=u*DDy'+by; 
    [wx,wy]=shrink2(temp1,temp2,1./sigma1);
    %bproblem
    bx=bx+DDx*u-wx;
    by=by+u*DDy'-wy;  
    
    % w2 problem
    temp3 = DDx*v + z*DDy' + b2;
    w2 = sign(temp3).*max(abs(temp3) - 1./sigma2,0);

    %bproblem
    b2=temp3-w2;  
    

    
    u = newu; v = newv; z = newz;
end

disp(['All iterations: ' num2str(i)]);

% linear streach, i.e. (3.8) in our paper
%min_u=min(min(u));
%max_u=max(max(u));
%u=(u-min_u)/(max_u-min_u);

%--------------------------------------------------------------------------
function [H1_l,H1_r] = first_order_differences(m,n)
%--------------------------------------------------------------------------
%derivatives at the right and lower boundary are set to zero (mirroring)

d1 = [-ones(m-1,1); 0];
d2 = ones(m,1);
H1_l = spdiags([d1,d2],[0,1],m,m);

d1 = [-ones(n-1,1); 0];
d2 = ones(n,1);
H1_r = spdiags([d1,d2],[0,1],n,n);



