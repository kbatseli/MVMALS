function yhat=sim_volterraTT(u,TT)
% yhat=sim_volterraTT(u,TT)
% --------------------------
% Simulates a truncated MISO Volterra series in the Tensor Train (TT) format for
% given inputs u(:,1),u(:,2),....
%
% yhat      =   vector, column vector of simulated output signal,
%
% u         =   matrix, u(:,k) contains the kth input series,
%
% TT        =   cell, TT{i} contains the ith TT core of the Volterra
%               tensor.
%
% Reference
% ---------
%
% 06-2016, Kim Batselier
p=size(u,2);                    % number of inputs
yhat=zeros(length(u),1);
d=length(TT);
r=zeros(1,d+1);
for i=1:d
    [r0,M,r1]=size(TT{i});
    M=(M-1)/p;
    r(i:i+1)=[r0 r1];
end
uextended=[zeros(M-1,p);u];     % append zeros to inputs      
for j=M:length(u)+M-1
    uj=zeros(p*M+1,1);
    uj(1)=1;
    for k=1:p
        uj(2+(k-1)*M:2+k*M-1)=uextended(j:-1:j-M+1,k);                
    end
    temp=uj'*reshape(TT{1},[p*M+1,r(2)]);    
    for k=2:d
        temp=temp*reshape(uj'*reshape(permute(TT{k},[2 3 1]),[p*M+1,r(k+1)*r(k)]),[r(k+1),r(k)])';
    end
    yhat(j-M+1)=temp;
end
end
