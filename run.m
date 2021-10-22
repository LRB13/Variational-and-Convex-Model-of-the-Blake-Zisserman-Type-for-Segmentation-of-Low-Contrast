%%
load('ultrasound.mat')
%% add noise if needed

Img_clean = im;
sigma = 10;
im = imnoise(im,'gaussian',0,(sigma/255)^2);

%%
%ultrasound
lambda1 = 15; mu1 = 1; sigma1 = 2; xi1=1000;
lambda2 = 20; mu2 = 1; sigma2 = 2; xi2=1000;
lambda3 = 20; mu3 = 1; sigma3 = 2; xi3=1000;


[u,v,z]=game_threeplayer(im,lambda1,mu1,lambda2,mu2,lambda3,mu3,sigma1,sigma2,sigma3,xi1,xi2,xi3);
%[u,v,z]=game_twoplayer(im,lambda1,mu1,lambda2,mu2,sigma1,sigma2,xi1,xi2);


figure;
subplot(1,3,1); imagesc(u); title("u, \lambda_1 = " + lambda1 + ", \mu_1 = " + mu1 + ", \xi_1 = " + xi1); colormap gray; axis image; colorbar;
subplot(1,3,2); imagesc(v); title("u_x, \lambda_2 = " + lambda2 + ", \mu_2 = " + mu2 + ", \xi_2 = " + xi2); colormap gray; axis image; colorbar;
subplot(1,3,3); imagesc(z); title("u_y, \lambda_3 = " + lambda3 + ", \mu_3 = " + mu3 + ", \xi_3 = " + xi3); colormap gray; axis image; colorbar

%s = sprintf('print -djpeg temp.jpg;'); eval(s)
grad = sqrt(v.^2 + z.^2);
%% k-means
k = 2;

%th=ThdKmeans(u,k); %define using kmeans
th = 0.48; %manual define

seg=SegResultShow(im,u,th,k); %this produces figures to show the segmentation result
bin_seg = seg>th; %for two phase segmentation

