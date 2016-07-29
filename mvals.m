function [TN,e]=mvals(y,u,M,r,varargin)
% [TN,e]=mvals(y,u,M,r,THRESHOLD)
% -------------------------------
% MIMO Volterra Alternating Linear Scheme (MVALS) algorithm for 
% solving the MIMO Volterra system identification problem in the Tensor
% Network format.
%
% TN        =   cell, TN{i} contains the ith TN core of the MIMO Volterra
%               tensor,
%
% e         =   vector, e(i) contains the relative residual 
% 				||y-yhat||_2/||y||_2 at iteration i,
%
% y         =   matrix, y(:,k) contains the kth output,
%
% u         =   matrix, u(:,k) contains the kth input,
%
% M         =   scalar, memory of each of the Volterra kernels,
%
% r         =   vector, contains the TT ranks r_1 up to r_{d-1}, since
%               r_0=r_d=1.
%
% THRESHOLD =   scalar, optional threshold on RMS error to stop iterations.
%               Default=1e-4.    
%
% Reference
% ---------
%
% 06-07-2016, Kim Batselier

p=size(u,2);                    % number of inputs
uextended=[zeros(M-1,p);u];     % append zeros to inputs  
[N,l]=size(y);
y=reshape(y',[N*l,1]);
d=length(r)+1;                  % degree of truncated Volterra series
r=[l r(:)' 1];                  % append extremal TN ranks
MAXITR=100;
if ~isempty(varargin)
    THRESHOLD=varargin{1};
else
    THRESHOLD=1e-4;
end 

% initialize right-orthonormal cores with prescribed TN ranks
TN=cell(1,d);
TN{1}=rand(r(1),p*M+1,r(2));
TN{1}=TN{1}./norm(TN{1}(:));
for i=2:d
    TN{i}=reshape(orth(rand((p*M+1)*r(i+1),r(i))),[r(i),(p*M+1),r(i+1)]);    
end

yhat=sim_volterraTN(u,TN);
yhat=reshape(yhat',[N*l,1]);
e(1)=norm(y(l*M+1:end)-yhat(l*M+1:end))/norm(y(l*M+1:end));
%e(1)=sqrt(norm(y-yhat)^2/N);    % RMS metric for residual

itr=1;                          % counts number of iterations
ltr=1;                          % flag that checks whether we sweep left to right
sweepindex=1;                   % index that indicates which TT core will be updated

while itr<2 || ((e(itr) < e(itr-1)) && (itr < MAXITR) && e(itr) > THRESHOLD)
    updateTT;
    updatesweep;
    % only check residual after 1 half sweep
    if (sweepindex==d) || (sweepindex==1) % half a sweep
        itr=itr+1;
        yhat=sim_volterraTN(u,TN);
        yhat=reshape(yhat',[N*l,1]);
        e(itr)=norm(y(l*M+1:end)-yhat(l*M+1:end))/norm(y(l*M+1:end));
%        e(itr)=sqrt(norm(y-yhat)^2/N);    % RMS metric for residual
    end    
end  

    function updateTT
        % first construct the linear subsystem matrix
        if N < (p*M+1)*prod(r(sweepindex:sweepindex+1))
            warning(['Matrix of reduced system is underdetermined: ' num2str(N) ' < ' num2str((p*M+1)*prod(r([sweepindex,sweepindex+2])))])
        end
        A=zeros(N,(p*M+1)*prod(r(sweepindex:sweepindex+1)));
        for i=M:N+M-1
            ui=zeros(p*M+1,1);
            ui(1)=1;
            for j=1:M
                ui(2+(j-1)*p:2+j*p-1)=uextended(i-j+1,:)';                
            end        
            vk1=eye(l);     % initialize row vector v_{k-1}
            vk2=1;          % initialize column vector v_{k+1}
            for j=1:sweepindex-1
                vk1=vk1*reshape(reshape(permute(TN{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            for j=sweepindex+1:d
                vk2=vk2* reshape(reshape(permute(TN{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            A((i-M)*l+1:(i-M+1)*l,:)=mkron(vk2',ui',vk1);
        end
        g=pinv(A)*y;
        if ltr
            % left-to-right sweep, generate left orthogonal cores
            [Q,R]=qr(reshape(g,[r(sweepindex)*(p*M+1),r(sweepindex+1)])); 
            TN{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),p*M+1,r(sweepindex+1)]);
            TN{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TN{sweepindex+1},[r(sweepindex+1),(p*M+1)*r(sweepindex+2)]),[r(sweepindex+1),p*M+1,r(sweepindex+2)]);
        else
            % right-to-left sweep, generate right orthogonal cores
            [Q,R]=qr(reshape(g,[r(sweepindex),(p*M+1)*r(sweepindex+1)])'); 
            TN{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),p*M+1,r(sweepindex+1)]);
            TN{sweepindex-1}=reshape(reshape(TN{sweepindex-1},[r(sweepindex-1)*(p*M+1),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),p*M+1,r(sweepindex)]);
        end
    end


    function updatesweep
        if ltr
            sweepindex=sweepindex+1;
            if sweepindex== d                
                ltr=0;
            end
        else
            sweepindex=sweepindex-1;
            if sweepindex== 1                
                ltr=1;
            end
        end
    end
end
