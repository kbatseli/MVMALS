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
% TN        =   cell, TN{i} contains the ith TN core of the MIMO Volterra
%               tensor.
%
% Reference
% ---------
%
% 06-2016, Kim Batselier
p=size(u,2);                    % number of inputs
yhat=zeros(length(u),1);
d=length(TN);
r=zeros(1,d+1);
for i=1:d
    [r0,M,r1]=size(TN{i});
    M=(M-1)/p;
    r(i:i+1)=[r0 r1];
end
uextended=[zeros(M-1,p);u];     % append zeros to inputs      
for j=M:length(u)+M-1
    uj=zeros(p*M+1,1);
    uj(1)=1;
    for k=1:M
        uj(2+(k-1)*p:2+k*p-1)=uextended(j-k+1,:)';                
    end
    temp=reshape(uj'*reshape(permute(TN{1},[2 1 3]),[(p*M+1),r(1)*r(2)]),[r(1),r(2)]);    
    for k=2:d
        temp=temp*reshape(uj'*reshape(permute(TN{k},[2 3 1]),[p*M+1,r(k+1)*r(k)]),[r(k+1),r(k)])';
    end
    yhat((j-M)*r(1)+1:(j-M+1)*r(1))=temp;
end
yhat=reshape(yhat,[r(1),length(u)])';
end