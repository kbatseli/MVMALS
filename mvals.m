function [TT,e]=mvals(y,u,M,r,varargin)
% [TT,e]=mvals(y,u,M,r,THRESHOLD)
% -------------------------------
% MIMO Volterra Alternating Linear Scheme (MVALS) algorithm for 
% solving the MIMO Volterra system identification problem in the Tensor
% Train format. This function handles the MISO case, use it for each output
% separately.

%
% TT        =   tensor, Volterra tensor as a Tensor Train,
%
% e         =   vector, e(i) contains the relative residual during iteration i
%
% y         =   vector, contains the output series,
%
% u         =   vector, contains the input series,
%
% M         =   scalar, memory of each of the Volterra kernels,
%
% r         =   vector, contains the TT ranks r_1 up to r_{d-1}, since
%               r_0=r_d=1.
%
% THRESHOLD =   scalar, optional threshold on relative residual to stop 
%               iterations. Default=1e-9.    
%
% Reference
% ---------
%
% 06-2016, Kim Batselier

d=length(r)+1;  % degree of truncated Volterra series
r=[1 r(:)' 1];  % append extremal ones to TT ranks
MAXITR=100;
if ~isempty(varargin)
    THRESHOLD=varargin{1};
else
    THRESHOLD=1e-4;
end

% make sure y,u are column vectors
p=size(u,2);                    % number of inputs
y=y(:);                         % make sure y is a column vector
uextended=[zeros(M-1,p);u];     % append zeros to inputs      
N=length(y);                    % number of samples

% initialize right-orthonormal cores with prescribed TT ranks
TT=cell(1,d);
TT{1}=rand(r(1),p*M+1,r(2));
TT{1}=TT{1}./norm(TT{1}(:));
for i=2:d
    TT{i}=reshape(orth(rand((p*M+1)*r(i+1),r(i))),[r(i),(p*M+1),r(i+1)]);    
end

yhat=sim_volterraTT(u,TT);
e(1)=norm(y-yhat)/norm(y); 

itr=1;                          % counts number of iterations
ltr=1;                          % flag that checks whether we sweep left to right
sweepindex=1;                   % index that indicates which TT core will be updated

while itr<2 || ((e(itr) < e(itr-1)) && (itr < MAXITR) && e(itr) > THRESHOLD)
    updateTT;
    updatesweep;
    % only check residual after 1 half sweep
    if (sweepindex==d) || (sweepindex==1) % half a sweep
        itr=itr+1;
        yhat=sim_volterraTT(u,TT);
        e(itr)=norm(y-yhat)/norm(y);     
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
            for j=1:p
                ui(2+(j-1)*M:2+j*M-1)=uextended(i:-1:i-M+1,j);                
            end
            vk1=1;  % initialize row vector v_{k-1}
            vk2=1;  % initialize column vector v_{k+1}
            for j=1:sweepindex-1
                vk1=vk1*reshape(reshape(permute(TT{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            for j=sweepindex+1:d
                vk2=vk2* reshape(reshape(permute(TT{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            A(i-M+1,:)=mkron(vk2',ui',vk1);
        end
        g=pinv(A)*y;
        if ltr
            % left-to-right sweep, generate left orthogonal cores
            [Q,R]=qr(reshape(g,[r(sweepindex)*(p*M+1),r(sweepindex+1)])); 
            TT{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),p*M+1,r(sweepindex+1)]);
            TT{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TT{sweepindex+1},[r(sweepindex+1),(p*M+1)*r(sweepindex+2)]),[r(sweepindex+1),p*M+1,r(sweepindex+2)]);
        else
            % right-to-left sweep, generate right orthogonal cores
            [Q,R]=qr(reshape(g,[r(sweepindex),(p*M+1)*r(sweepindex+1)])'); 
            TT{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),p*M+1,r(sweepindex+1)]);
            TT{sweepindex-1}=reshape(reshape(TT{sweepindex-1},[r(sweepindex-1)*(p*M+1),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),p*M+1,r(sweepindex)]);
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