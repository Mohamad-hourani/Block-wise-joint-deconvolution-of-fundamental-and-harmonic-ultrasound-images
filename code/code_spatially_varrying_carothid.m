%=========================================================================%
%========   Joint deconvolution of fundanmental and harmonic      ========%
%========                ultrasound images                        ========%
%========           Code Author: Mohamad HOURANI                  ========%
%========          Version (Date): Feb,07 2020                    ========%
%========            Email: mohamad.hourani@irit.fr               ========%
%=========================================================================%
%------------------------      COPYRIGHT      ----------------------------%
% Copyright (2020): Mohamad HOURANI, Adrian Basarab, Denis Kouam\'{e}, and 
% Jean-Yves Tourneret;
% Permission to use, copy, modify, and distribute this software for any
% purpose without fee is hereby granted, provided that this entire notice
% is included in all copies of any software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or iplied
% warranty. Inparticular, the authors do not make any representation or
% warranty of any kind concerning the merchantability of this software or
% its fitness for any particular purpose.
%-------------------------------------------------------------------------%
%-----------------------      REFERENCE      -----------------------------%
% This set of MATLAB files contain an implementation of the algorithms
% described in the following paper:
%
% Mohamad Hourani, Adrian Basarab, Denis Kouam\'{e}, Jean-Yves Tourneret, 
% Guilia Matrone, Alessandro Ramalli,"On Ultrasound Image Deconvolution 
% using Fundamental and Harmonic Data"
%-------------------------------------------------------------------------%


%% add paths
addpath ./../data;
addpath ./../utils;
addpath ./../psf_est
%% Data 
load sgnFIG5;
rf=rf;
%% Attenuation (exponential)
[SF,SH,r,z]=spectrum_of_cells_gliss(rf./max(rf(:)),110,1,fs);
W=z;
Winv=1./z;
clear sgnFIG5;  
load hsgnFIG5;
load fsgnFIG5;
rff=rff;
rff2=rff2;
y1 = rff /max(rff(:));
y2 = rff2 /max(rff2(:));
load('h_est.mat','h_estf1','h_esth1','h_estf2','h_esth2','h_estf3','h_esth3','h_estf4','h_esth4',...
    'h_estf5','h_esth5','h_estf6','h_esth6','h_estf7','h_esth7','h_estf8','h_esth8')
%% -----------------------------------------------PSF1------------------------------------------------- 
% [h_estf1,psf1,phf1]=psf_est(y1(450:585,174:200),0.86);                                         % Estimation of the fundamental PSF
% [h_esth1,psh1,phh1]=psf_est(y2(450:585,174:200),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf1,1)/2);
H1f=h_estf1(n-40+1:n+40+1,n-10+1:n+10);
H1h=h_esth1(n-40+1:n+40+1,n-10+1:n+10+1);
H1f=H1f/max(H1f(:));
H1h=H1h/max(H1h(:));
H1=H1f;
H2=H1h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,
subplot(121),imagesc(xAxis,zAxis,rf2bmode(y1(1:3500,:),0.01));colormap(gray)%;title('Fundamental ultrasound image')
subplot(122),imagesc(xAxis,zAxis,rf2bmode(y2(1:3500,:),0.01));colormap(gray)%;title('Harmonic ultrasound image Image (with attenuation)')
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x1,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint1= x1./max(x1(:));
l=rf2bmode(x_joint1,1);
lx1=20*log10(abs(hilbert(x_joint1)));
%% PSF2
% [h_estf2,psf2,phf2]=psf_est(y1(1025:1180,175:200),0.86);                                         % Estimation of the fundamental PSF
% [h_esth2,psh2,phh2]=psf_est(y2(1025:1180,175:200),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf2,1)/2);
H2f=h_estf2(n-40+1:n+40+1,n-10+1:n+10);
H2h=h_esth2(n-40+1:n+40+1,n-10+1:n+10+1);
H2f=H2f/max(H2f(:));
H2h=H2h/max(H2h(:));
H1=H2f;
H2=H2h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x2,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint2= x2./max(x2(:));
l=rf2bmode(x_joint2,1);
lx2=20*log10(abs(hilbert(x_joint2)));
%% PSF3
% [h_estf3,psf3,phf3]=psf_est(y1(1555:1705,177:200),0.86);                                         % Estimation of the fundamental PSF
% [h_esth3,psh3,phh3]=psf_est(y2(1555:1705,177:200),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf3,1)/2);
H3f=h_estf3(n-40+1:n+40+1,n-10+1:n+10);
H3h=h_esth3(n-40+1:n+40+1,n-10+1:n+10+1);
H3f=H3f/max(H3f(:));
H3h=H3h/max(H3h(:));
H1=H3f;
H2=H3h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x3,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint3= x3./max(x3(:));
l=rf2bmode(x_joint3,1);
lx3=20*log10(abs(hilbert(x_joint3)));
%% PSF4
% [h_estf4,psf4,phf4]=psf_est(y1(2070:2225,179:205),0.86);                                         % Estimation of the fundamental PSF
% [h_esth4,psh4,phh4]=psf_est(y2(2070:2225,179:205),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf4,1)/2);
H4f=h_estf4(n-40+1:n+40+1,n-10+1:n+10);
H4h=h_esth4(n-40+1:n+40+1,n-10+1:n+10+1);
H4f=H4f/max(H4f(:));
H4h=H4h/max(H4h(:));
H1=H4f;
H2=H4h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x4,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint4= x4./max(x4(:));
l=rf2bmode(x_joint4,1);
lx4=20*log10(abs(hilbert(x_joint4)));
%% PSF5
% [h_estf5,psf5,phf5]=psf_est(y1(2603:2747,179:202),0.86);                                         % Estimation of the fundamental PSF
% [h_esth5,psh5,phh5]=psf_est(y2(2603:2747,179:202),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf5,1)/2);
H5f=h_estf5(n-40+1:n+40+1,n-10+1:n+10);
H5h=h_esth5(n-40+1:n+40+1,n-10+1:n+10+1);
H5f=H5f/max(H5f(:));
H5h=H5h/max(H5h(:));
H1=H5f;
H2=H5h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x5,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint5= x5./max(x5(:));
l=rf2bmode(x_joint5,1);
lx5=20*log10(abs(hilbert(x_joint5)));
%% PSF6
% [h_estf6,psf6,phf6]=psf_est(y1(3101:3239,180:203),0.86);                                         % Estimation of the fundamental PSF
% [h_esth6,psh6,phh6]=psf_est(y2(3101:3239,180:203),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf6,1)/2);
H6f=h_estf6(n-40+1:n+40+1,n-10+1:n+10);
H6h=h_esth6(n-40+1:n+40+1,n-10+1:n+10+1);
H6f=H6f/max(H6f(:));
H6h=H6h/max(H6h(:));
H1=H6f;
H2=H6h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x6,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint6= x6./max(x6(:));
l=rf2bmode(x_joint6,1);
lx6=20*log10(abs(hilbert(x_joint6)));
%% PSF7
% [h_estf7,psf7,phf7]=psf_est(y1(3594:3727,181:204),0.86);                                         % Estimation of the fundamental PSF
% [h_esth7,psh7,phh7]=psf_est(y2(3594:3727,181:204),0.86);                                         % Estimation of the harmonic PSF[h_est1,ps,ph]=psf_est(y1(100:1100,100:200),0.86);                                         % Estimation of the fundamental PSF
n=round(size(h_estf7,1)/2);
H7f=h_estf7(n-40+1:n+40+1,n-10+1:n+10);
H7h=h_esth7(n-40+1:n+40+1,n-10+1:n+10+1);
H7f=H7f/max(H7f(:));
H7h=H7h/max(H7h(:));
H1=H7f;
H2=H7h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x7,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint7= x7./max(x7(:));
l=rf2bmode(x_joint7,1);
lx7=20*log10(abs(hilbert(x_joint7)));
%% PSF8 
% [h_estf8,psf8,phf8]=psf_est(y1(4100:4218,182:204),0.86);                                         % Estimation of the fundamental PSF
% [h_esth8,psh8,phh8]=psf_est(y2(4100:4218,182:204),0.86);                                         % Estimation of the harmonic PSF
n=round(size(h_estf8,1)/2);
H8f=h_estf8(n-40+1:n+40+1,n-10+1:n+10);
H8h=h_esth8(n-40+1:n+40+1,n-10+1:n+10+1);
H8f=H8f/max(H8f(:));
H8h=H8h/max(H8h(:));
H1=H8f;
H2=H8h;
% %% Preparation of inputs 
N = numel(y1);
BSNRdb1 = 40;                                                              % Blur signal to noise ratio
BSNRdb2 = 40;                                                              % Blur signal to noise ratio
sigma1 = norm(y1-mean(mean(y1)),'fro')/sqrt(N*10^(BSNRdb1/10));            % Standard deviation of the WGN in fundamental image
sigma2 = norm(y2-mean(mean(y2)),'fro')/sqrt(N*10^(BSNRdb2/10));            % Standard deviation of the WGN in fundamental image
%% PSF estimation and operators
[FB1,FBC1,F2B1,Bpad1] = Operator_PSF(H1,y1);                               % Operator of the fundamental PSF
[FB2,FBC2,F2B2,Bpad2] = Operator_PSF(H2,y2);                               % Operator of the harmonic PSF
% Normalized data
%% Presenting inputs
figure,subplot(121),imagesc(H1),title('Fundamental PSF')
subplot(122),imagesc(H2),title('Harmonic PSF')
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----Joint Solution----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2*||y1-H1x||^2 +1/2*||y2-WH2x||^2 + mu||x||_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of ADMM 
stoppingrule = 1;                                                          % stopping rule (1- difference between two objectives,  2- distances(norm value)
tolA = -inf;                                                               % Tolerance for stopping criterion
sig2n = sigma1^2;                                                          % Variance of noise
p =1;                                                                      % lp normin the prior
beta =5;                                                                   % AL Hyperparameter
taup =5e-2;                                                                    
mu  = taup*sig2n; %mu=[]                                                   % lp norm hyperparameter
nt= numel(mu);                                                             % number of the set of lp norm hyperparameter
maxiter = 50;                                                              % maximum number of iteration                           
objective = zeros(nt,maxiter);                                             % objective's initialization
times = zeros(1,maxiter);                                                  % time of ADMM iterations
% Inputs
X1 = y1;
X2 = y2;
U = y1;
D = 0*y1;
%% ADMM Iterations
tic
[x8,U,D,criterion]=joint_deconv_ADMM(FB1,FB2,FBC1,FBC2,F2B1,F2B2,W,Winv,X1,X2,U,D,beta,mu,p,tolA,objective,times,stoppingrule,maxiter);
toc
%figure,plot(criterion);
x_joint8= x8./max(x8(:));
l=rf2bmode(x_joint8,1);
lx8=20*log10(abs(hilbert(x_joint8)));
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zaxis=linspace(0,1,size(zAxis,2));
a=asin(0.5)/0.1248;
y1=1+sin(a*zaxis+pi);
y1(561:end)=0;
y2=-4.0064*zaxis+0.5;
y2(560:end)=0;
w1=y1+circshift(y2,560);
figure,plot(w1)
y11=sin(a*zaxis);
y11(560:end)=0;
y3=zeros(size(zaxis));
y3(560:2*560)=0.5;
w2=y11+y3+circshift(y2,560*2)
figure,plot(w2);
y4=4.0064*zaxis;
y4(560:end)=0;
w3=circshift(y4,560)+circshift(y3,560)+circshift(y2,3*560);
figure,plot(w3)
w4=circshift(w3,560);
w5=circshift(w4,560);
w6=circshift(w5,560);
w7bis=circshift(w6,560);
w7bis(3921:end)=0;
w7=flip(w2);
w7=circshift(w7,1)
w8=flip(w1);
w8=circshift(w8,-1);
figure,plot(w1),hold on, plot(w2), hold on, plot(w3), hold on, plot(w4), hold on, plot(w5), hold on, plot(w6), hold on, plot(w7), hold on, plot(w8)
W1=repmat(w1,size(rf,2),1);
W2=repmat(w2,size(rf,2),1);
W3=repmat(w3,size(rf,2),1);
W4=repmat(w4,size(rf,2),1);
W5=repmat(w5,size(rf,2),1);
W6=repmat(w6,size(rf,2),1);
W7=repmat(w7,size(rf,2),1);
W8=repmat(w8,size(rf,2),1);
x_joint_spatially = W1'.*x_joint1 + W2'.*x_joint2 + W3'.*x_joint3 + W4'.*x_joint4...
    + W5'.*x_joint5 + W6'.*x_joint6 + W7'.*x_joint7 + W8'.*x_joint8;
load('x_joint_carotid.mat');
figure,subplot(121),imagesc(xAxis,zAxis,rf2bmode(x_joint(1:3500,:),0.01)),colormap('gray')
subplot(122),imagesc(xAxis,zAxis,rf2bmode(x_joint_spatially(1:3500,:),0.01)),colormap('gray')
figure,subplot(121),imagesc(xAxis,zAxis,rf2bmode(x_joint(100:1000,:),0.01)),colormap('gray')
subplot(122),imagesc(xAxis,zAxis,rf2bmode(x_joint_spatially(100:1000,:),0.01)),colormap('gray')
% %% Profile horizantal
% % PSF1
%  figure,plot(rf2bmode(x_joint(514,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(514,:),0.1)),legend
% %PSF2
%  figure,plot(rf2bmode(x_joint(1085,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(1085,:),0.1)),legend
% %PSF3
%  figure,plot(rf2bmode(x_joint(1622,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(1622,:),0.1)),legend
% %PSF4
% figure,plot(rf2bmode(x_joint(2137,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(2137,:),0.1)),legend
% %PSF5
% figure,plot(rf2bmode(x_joint(2674,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(2674,:),0.1)),legend
% %PSF6
% figure,plot(rf2bmode(x_joint(3157,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(3157,:),0.1)),legend
% %PSF7
% figure,plot(rf2bmode(x_joint(3664,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(3664,:),0.1)),legend
% %PSF8
% figure,plot(rf2bmode(x_joint(4147,:),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(4147,:),0.1)),legend
% %% Profile vertical
% % PSF1
%  figure,plot(rf2bmode(x_joint(450:700,187),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(450:700,187),0.1)),legend
% %PSF2
%  figure,plot(rf2bmode(x_joint(800:1200,188),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(800:1200,188),0.1)),legend
% %PSF3
%  figure,plot(rf2bmode(x_joint(1400:1800,192),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(1400:1800,192),0.1)),legend
% %PSF4
% figure,plot(rf2bmode(x_joint(1900:2300,193),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(1900:2300,193),0.1)),legend
% %PSF5
% figure,plot(rf2bmode(x_joint(2400:2800,191),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(2400:2800,191),0.1)),legend
% %PSF6
% figure,plot(rf2bmode(x_joint(2900:3300,192),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(2900:3300,192),0.1)),legend
% %PSF7
% figure,plot(rf2bmode(x_joint(3400:3800,194),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(3400:3800,194),0.1)),legend
% %PSF8
% figure,plot(rf2bmode(x_joint(3900:4300,195),0.1))
% , hold on, plot(rf2bmode(x_joint_spatially(3900:4300,195),0.1)),legend
