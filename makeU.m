function U=makeU(u,M,d)
% U=makeU(u,M,d)
% --------------
% Generates the N x (pM+1)^d matrix U such that y=U*vec(V) for a p-input
% Volterra system, where vec(V) is the vectorization of the MIMO Volterra
% tensor and N=size(u,1) is the number of measured samples. 
%
% U 	=   matrix, N x (pM+1)^d matrix U such that y=U*vec(V),
%
% u     =   matrix, u(:,k) contains the kth input series,
%
% M		=	scalar, memory length of the Volterra system,
%
% d 	=	scalar, maximal degree of the truncated Volterra series. 	
%
% Reference
% ---------
%
% 07/2016, Kim Batselier


[N,p]=size(u);
U=zeros(N,(p*M+1)^d);
u=[zeros(M-1,p);u];
for i=M:N+M-1            
	temp=ones(1,p*M+1);
	for j=1:p
		temp(2+(j-1)*M:2+j*M-1)=u(i:-1:i-M+1,j)';
	end
	U(i-M+1,:)=mkron(temp,d);
end

end
