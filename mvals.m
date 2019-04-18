function [TN,e]=mvals(y,u,M,r,varargin)
% [TN,e]=mvals(y,u,M,r,THRESHOLD)
% -------------------------------
% MIMO Volterra Alternating Linear Scheme (MVALS) algorithm for 
% solving the MIMO Volterra system identification problem in the Tensor
% Network format.
%
% TN        =   Tensor Network, TN.core is a cell containing the TN-cores,
%               TN.n is a matrix where TN.n(i,:) are the dimensions of the
%               ith TN-core,
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
% 06/07/11 - 2016, Kim Batselier

p=size(u,2);                    % number of inputs
uextended=[zeros(M-1,p);u];     % append zeros to inputs  
[N,l]=size(y);
y=reshape(y',[N*l,1]);
d=length(r)+1;                  % degree of truncated Volterra series
r=[l r(:)' 1];                  % append extremal TN ranks
MAXITR=100;
n=p*M+1;
if ~isempty(varargin)
    THRESHOLD=varargin{1};
else
    THRESHOLD=1e-4;
end

% construct N x n matrix U
U=zeros(N,n);
u=[zeros(M-1,p);u];
for i=M:N+M-1            
	temp=ones(1,n);
    for j=1:M
        temp(2+(j-1)*p:2+j*p-1)=u(i-j+1,:);                
    end   
    U(i-M+1,:)=temp;
end
u=u(M:end,:);
Vp=cell(1,d);
Vm=cell(1,d);
if l==1
    Vm{1}=ones(N,1);
else
    Vm{1}=eye(l);
end
Vp{d}=ones(N,1);

% initialize right-orthonormal cores with prescribed TN ranks
TN.core=cell(1,d);
TN.core{1}=rand(r(1),n,r(2));
TN.core{1}=TN.core{1}./norm(TN.core{1}(:));
TN.n(1,:)=[1 l n r(2)];
for i=d:-1:2
	TN.n(i,:)=[r(i) 1 n r(i+1)];
    TN.core{i}=permute(reshape(orth(rand((n)*r(i+1),r(i))),[r(i+1),(n),r(i)]),[3,2,1]);    
    Vp{i-1}=dotkron(Vp{i},U)*reshape(permute(TN.core{i},[3 2 1]),[r(i+1)*n,r(i)]); % N x r_{i-1}    
end

yhat=sim_volterraTN(u,TN);
yhat=reshape(yhat',[N*l,1]);
e(1)=norm(y(l*M+1:end)-yhat(l*M+1:end))/norm(y(l*M+1:end));

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
    end    
end  

    function updateTT
%         % first construct the linear subsystem matrix
        if l==1
            A=dotkron(Vm{sweepindex},U,Vp{sweepindex});
        elseif sweepindex == 1
            A=kron(dotkron(U,Vp{sweepindex}),Vm{sweepindex});
        else
            A=dotkron(Vm{sweepindex},U,Vp{sweepindex});
            A=reshape(A,[N,l,r(sweepindex)*n*r(sweepindex+1)]);
            A=permute(A,[2 1 3]);
            A=reshape(A,[N*l,r(sweepindex)*n*r(sweepindex+1)]);
        end 
        g=pinv(A)*y;
        if ltr
            % left-to-right sweep, generate left orthogonal cores and update vk1
            [Q,R]=qr(reshape(g,[r(sweepindex)*(n),r(sweepindex+1)])); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),n,r(sweepindex+1)]);
            TN.core{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TN.core{sweepindex+1},[r(sweepindex+1),(n)*r(sweepindex+2)]),[r(sweepindex+1),n,r(sweepindex+2)]);
            if l==1
                Vm{sweepindex+1}=dotkron(Vm{sweepindex},U)*reshape(TN.core{sweepindex},[r(sweepindex)*n,r(sweepindex+1)]); % N x r_{i}
            elseif sweepindex==1
                Vm{sweepindex+1}=U*reshape(permute(TN.core{sweepindex},[2 1 3]),[n,r(sweepindex)*r(sweepindex+1)]); %N x r_{i-1}r_i
                Vm{sweepindex+1}=reshape(Vm{sweepindex+1},[N,r(sweepindex)*r(sweepindex+1)]);                
            else
                Vm{sweepindex+1}=reshape(dotkron(Vm{sweepindex},U),[N*l,r(sweepindex)*n])*reshape(TN.core{sweepindex},[r(sweepindex)*n,r(sweepindex+1)]);
                Vm{sweepindex+1}=reshape(Vm{sweepindex+1},[N,l*r(sweepindex+1)]);
            end
        else
            % right-to-left sweep, generate right orthogonal cores and update vk2
            [Q,R]=qr(reshape(g,[r(sweepindex),(n)*r(sweepindex+1)])'); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),n,r(sweepindex+1)]);
            TN.core{sweepindex-1}=reshape(reshape(TN.core{sweepindex-1},[r(sweepindex-1)*(n),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),n,r(sweepindex)]);
            Vp{sweepindex-1}=dotkron(Vp{sweepindex},U)*reshape(permute(TN.core{sweepindex},[3 2 1]),[r(sweepindex+1)*n,r(sweepindex)]); % N x r_{i-1}    

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
