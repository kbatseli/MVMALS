%% SISO decaying exponentials
clear all
load decexp

% determine maximal TT-ranks for 700 samples with MALS
tic;
[TT,e]=mvmals(y(1:700),u(1:700),M,d);
toc
% verify on validation data
yhat=sim_volterraTN(u,TT);
plot(y(701:720),'-go');grid on,hold on,plot(yhat(701:720))

% compare runtime with ALS for same number of samples and TT-ranks
r=TNranks(TT);    
tic;
[TT2,e2]=mvals(y(1:700,:),u(1:700,:),M,r(2:end-1));
toc
yhat2=sim_volterraTN(u,TT2);
plot(yhat2(701:720))

%% MIMO example: 
%  3-input-2-output Volterra system with M=3,d=3

clear all
load simdata       

% use ALS for the first 1100 samples and fixed TN-ranks
[TN,e]=mvals(y(1:floor(1.1*(p*M+1)^d),:),u(1:floor(1.1*(p*M+1)^d),:),M,[11 10]);
% check relative residual 
e(end)
% get TN-ranks
r=TNranks(TN);  

% use MALS to automatically determine the TN-ranks
[TN2,e2]=mvmals(y(1:floor(1.1*(p*M+1)^d),:),u(1:floor(1.1*(p*M+1)^d),:),M,d);
e2(end)     % similar relative residual error but larger TN-ranks
% get TN-ranks
r2=TNranks(TN2);  
[r;r2]