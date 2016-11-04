function yhat=sim_volterraTN(u,TN)
% yhat=sim_volterraTN(u,TN)
% --------------------------
% Simulates a truncated MIMO Volterra series in the Tensor Network (TN) format for
% given inputs u(:,1),u(:,2),....
%
% yhat      =   matrix, y(:,k) contains the kth simulated output,
%
% u         =   matrix, u(:,k) contains the kth input,
%
% TN        =   Tensor Network.
%
% Reference
% ---------
%
% 06/11-2016, Kim Batselier

% generate N x (pM+1) matrix U with input samples
[N,p]=size(u);
d=length(TN.core);
n=TN.n(1,3);
l=TN.n(1,2);
M=(n-1)/p;
U=zeros(N,n);
u=[zeros(M-1,p);u];
for i=M:N+M-1            
	temp=ones(1,n);
    for j=1:M
        temp(2+(j-1)*p:2+j*p-1)=u(i-j+1,:);                
    end   
    U(i-M+1,:)=temp;
end

yhat=U*reshape(permute(reshape(TN.core{1},TN.n(1,:)),[3 1 2 4]),[n,prod(TN.n(1,[1 2 4]))]);
for i=2:d
    temp=dotkron(yhat,U);
    yhat=reshape(temp,[N*l,TN.n(i,1)*n])*reshape(TN.core{i},[prod(TN.n(i,1:3)),TN.n(i,end)]);
    yhat=reshape(yhat,[N,l*TN.n(i,end)]);
end
yhat=reshape(yhat,[N,l]);

end