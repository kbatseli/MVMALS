% Extended Complex Kalman Filter as described in
% Nonlinear system modeling and identification using Volterra-PARAFAC models
% by G. Favier, A.Y. Kibangou and T. Bouilloc

% first setup SISO system
clear all;
d=5;				% degree of SISO Volterra system
M=5;				% memory of SISO Volterra system, according to Volterra-PARAFAC definition
N=5000;				% total number of samples
Nsysid=4700;		% number of samples used for identification
u=randn(N,1);		% input
y=zeros(N,1);		% output
sigma_e=1e-6;		% measurement noise standard deviation
alpha=1e1;			% variance initial guess unknown parameters
lambda=1e3;			% forgetting factor

% create random factor matrices
R=M;
A{1}=randn(1,M);
for i=2:d
	A{i}=randn(R,M);
end

% simulate the system		
for i=M+1:N
	uk=u(i-1:-1:i-M);
	for j=1:d,y(i)=y(i)+sum((A{j}*uk).^j);end
end
y=y+sigma_e*randn(N,1);

% initialization
theta=randn(1+M+(d-1)*M*R,1);
K=length(theta);
P=alpha*eye(K);

% running the filter
for i=M+1:Nsysid
	uk=u(i-1:-1:i-M);
	g=zeros(K,1);
	psi=zeros(K,1);
	g(1:M+1)=[1;uk];
	psi(1:M+1)=[1;uk];
	for j=2:d
		% extract j-th factor matrix from theta
		Ahat=reshape(theta(M+2+(j-2)*M*R:M+1+(j-1)*M*R),[R,M]);
		v=(Ahat*uk).^(j-1);
		phi=j*kron(v,uk);
		g(M+2+(j-2)*M*R:M+1+(j-1)*M*R)=phi;
		psi(M+2+(j-2)*M*R:M+1+(j-1)*M*R)=phi/j;
	end
	e=y(i)-psi'*theta;
	sigma=g'*P*g+sigma_e^2;
	G=P*g/sigma;
	theta=theta+G*e;
	P=(P-G*G'*sigma)/lambda;
end

% extract the estimated factor matrices and simulate the system
clear Ahat
Ahat{1}=theta(2:M+1)';
for j=2:d
	Ahat{j}=reshape(theta(M+2+(j-2)*M*R:M+1+(j-1)*M*R),[R,M]);
end
yhat=zeros(N,1);
for i=M+1:N
	uk=u(i-1:-1:i-M);
	for j=1:d,
		yhat(i)=yhat(i)+theta(1)+sum((Ahat{j}*uk).^j);
	end
end

[TN,e]=mvals(y(1:700,:),u(1:700,:),M+1,(M+1)*ones(1,d-1));
yhat2=sim_volterraTN(u,TN);
norm(y(Nsysid+1:end)-yhat2(Nsysid+1:end))/norm(y(Nsysid+1:end))
