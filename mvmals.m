function [TT,e]=mvmals(y,u,M,d,varargin)
% [TT,e]=mvmals(y,u,M,d,THRESHOLD)
% -------------------------------
% MIMO Volterra Modified Alternating Linear Scheme (MVMALS) algorithm for 
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
% u         =   matrix, u(:,k) contains the kth input series,
%
% M         =   scalar, memory of each of the Volterra kernels,
%
% d         =   scalar, degree of truncated Volterra series, minimal d=2.
%
% THRESHOLD =   scalar, optional threshold on relative residual to stop
%               iterations. Default=1e-4.    
%
% Reference
% ---------
%
% 06-2016, Kim Batselier

MAXITR=10;
if ~isempty(varargin)
    THRESHOLD=varargin{1};
else
    THRESHOLD=1e-4;
end

p=size(u,2);                    % number of inputs
y=y(:);                         % make sure y is a column vector
uextended=[zeros(M-1,p);u];     % append zeros to inputs      
N=length(y);                    % number of samples
r=ones(1,d+1);                  % initial TT ranks to all 1's

% initialize random cores
TT=cell(1,d);
for i=1:d
    TT{i}=rand(r(i),p*M+1,r(i+1));
    TT{i}=TT{i}/norm(TT{i}(:));
end
yhat=sim_volterraTT(u,TT);
e(1)=norm(y-yhat)/norm(y); 

itr=1;                          % counts number of half-sweeps
ltr=1;                          % flag that checks whether we sweep left to right or right to left
sweepindex=1;                   % index that indicates which TT core will be updated

while (itr <2) || ( (e(itr) < e(itr-1)) && (itr < MAXITR) && e(itr) > THRESHOLD)
    updateTT;
    if (sweepindex==d-1 && ltr) || (sweepindex==1 && ~ltr)          % check whether half a sweep passed
        itr=itr+1;
        yhat=sim_volterraTT(u,TT);
        e(itr)=norm(y-yhat)/norm(y);     
    end     
    updatesweep;    
end  

    function updateTT
        % first construct the linear subsystem matrix
        if N < (p*M+1)^2*prod(r([sweepindex,sweepindex+2]))
            warning(['Matrix of reduced system is underdetermined: ' num2str(N) ' < ' num2str((p*M+1)^2*prod(r([sweepindex,sweepindex+2])))])
        end
        A=zeros(N,(p*M+1)^2*prod(r([sweepindex,sweepindex+2])));
        for i=M:N+M-1
            ui=zeros(p*M+1,1);
            ui(1)=1;
            for j=1:p
                ui(2+(j-1)*M:2+j*M-1)=uextended(i:-1:i-M+1,j);                
            end
            vk1=1;  % initialize row vector v_{k-1}
            vk2=1;  % initialize column vector v_{k+2}
            for j=1:sweepindex-1
               vk1=vk1*reshape(reshape(permute(TT{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            for j=sweepindex+2:d
                vk2=vk2* reshape(reshape(permute(TT{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end        
            A(i-M+1,:)=mkron(vk2',mkron(ui',2),vk1);
        end
        g=pinv(A)*y;                                        % truncated SVD solution (pseudoinverse)
        g=reshape(g,[r(sweepindex)*(p*M+1),(p*M+1)*r(sweepindex+2)]);
        [Ut,St,Vt]=svd(g);
        tol=eps(St(1))*max(size(g));
        st=diag(St);
        rankg=sum(st > tol);
        r(sweepindex+1)=rankg;
        if ltr
            % left-to-right sweep, generate left orthogonal cores
            TT{sweepindex}=reshape(Ut(:,1:rankg),[r(sweepindex),(p*M+1),rankg]);
            TT{sweepindex+1}=reshape(St(1:rankg,1:rankg)*Vt(:,1:rankg)',[rankg,p*M+1,r(sweepindex+2)]);
        else
            % right-to-left sweep, generate right orthogonal cores
            TT{sweepindex}=reshape(Ut(:,1:rankg)*St(1:rankg,1:rankg),[r(sweepindex),(p*M+1),rankg]);
            TT{sweepindex+1}=reshape(Vt(:,1:rankg)',[rankg,p*M+1,r(sweepindex+2)]);
        end        
    end


    function updatesweep
        if ltr
            if sweepindex < d-1
                sweepindex=sweepindex+1;
            else % sweepindex has reached end of the TT
                ltr=0; 
            end
        else
            if sweepindex > 1
                sweepindex=sweepindex-1;
            else % sweepindex has reached beginning of the TT
                ltr=1;
            end
        end
    end
end