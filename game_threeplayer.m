function [u,v,z]=game_threeplayer(Img,lambda1,mu1,lambda2,mu2,lambda3,mu3,sigma1,sigma2,sigma3,xi1,xi2,xi3)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


%% switches
update_dynamically=1;
%%
% initialization
[n,m] = size(Img);
u=Img;wx=0*Img;wy=wx;bx=wx;by=wx;
v=zeros(n,m); w2x=v; w2y=v; b2x=v; b2y=v;
z=zeros(n,m); w3x=z; w3y=z; b3x=z; b3y=z;
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
uker3 = lambda3 + xi3 +(mu3+sigma3)*fft2(uker);


Imgx = DDx*Img; Imgy = Img*DDy';
[Imgy,Imgx] = gradient(Img);

%%
for i=1:300
    % compute the right hand side 
    u0=u;
    rhs1 = lambda1.*Img+sigma1*DDx'*(wx-bx)+sigma1*(wy-by)*DDy + xi1*DDx'*v + xi1*z*DDy;
    newu = ifft2(fft2(rhs1)./uker1);
    if update_dynamically == 1
    u = newu;
    end
    
    rhs2 = lambda2.*Imgx + xi2.*DDx*u +sigma2*DDx'*(w2x-b2x)+sigma2*(w2y-b2y)*DDy;
    newv = ifft2(fft2(rhs2)./uker2);
    if update_dynamically == 1
    v = newv;
    end
    
    
    rhs3 = lambda3.*Imgy + xi3*u*DDy' +sigma3*DDx'*(w3x-b3x)+sigma3*(w3y-b3y)*DDy;
    newz = ifft2(fft2(rhs3)./uker3);
    if update_dynamically == 1
        z = newz;
    end    
    
   
    
    
    
    err=norm(newu-u0,'fro')/norm(u,'fro');
    
    if mod(i,10)==0
        disp(['iterations: ' num2str(i) '!  ' 'error is:   ' num2str(err)]);
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
    
    % w problem

    temp1=DDx*v+b2x;temp2=v*DDy'+b2y; 
    [w2x,w2y]=shrink2(temp1,temp2,1./sigma2);
    %bproblem
    b2x=b2x+DDx*v-w2x;
    b2y=b2y+v*DDy'-w2y;  
    
    % w problem

    temp1=DDx*z+b3x;temp2=z*DDy'+b3y; 
    [w3x,w3y]=shrink2(temp1,temp2,1./sigma3);
    %bproblem
    b3x=b3x+DDx*z-w3x;
    b3y=b3y+z*DDy'-w3y;  
    
    
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



