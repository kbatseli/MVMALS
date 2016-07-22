function [TN,e]=mvmals(y,u,M,d,varargin)
% [TN,e]=mvmals(y,u,M,d,THRESHOLD)
% -------------------------------
% MIMO Volterra Modified Alternating Linear Scheme (MVMALS) algorithm for 
% solving the MIMO Volterra system identification problem in the Tensor
% Network format. Upper bounds on the TN ranks are limited by the total
% number of output samples.
%
% TN        =   cell, TN{i} contains the ith TN core of the MIMO Volterra
%               tensor,
%
% e         =   vector, e(i) contains the root mean square (RMS) error at
%               iteration i,
%
% y         =   matrix, y(:,k) contains the kth output,
%
% u         =   matrix, u(:,k) contains the kth input,
%
% M         =   scalar, memory of each of the Volterra kernels,
%
% d         =   scalar, degree of truncated Volterra series, minimal d=2.
%
% THRESHOLD =   scalar, optional threshold on RMS error to stop iterations.
%               Default=1e-4.    
%
% Reference
% ---------
%
% 06-2016, Kim Batselier

p=size(u,2);                    % number of inputs
uextended=[zeros(M-1,p);u];     % append zeros to inputs  
[N,l]=size(y);
y=reshape(y',[N*l,1]);
r=[l ones(1,d)];                % initialize TN ranks to l and all 1's

MAXITR=10;
if ~isempty(varargin)
    THRESHOLD=varargin{1};
else
    THRESHOLD=1e-4;
end

% initialize random cores
TN=cell(1,d);
for i=1:d
    TN{i}=rand(r(i),p*M+1,r(i+1));
    TN{i}=TN{i}/norm(TN{i}(:));
end
yhat=sim_volterraTN(u,TN);
yhat=reshape(yhat',[N*l,1]);
e(1)=sqrt(norm(y-yhat)^2/N);    % RMS metric for residual

itr=1;                          % counts number of half-sweeps
ltr=1;                          % flag that checks whether we sweep left to right or right to left
sweepindex=1;                   % index that indicates which TT core will be updated

while (itr <2) || ( (e(itr) < e(itr-1)) && (itr < MAXITR) && e(itr) > THRESHOLD)
    updateTT;
    if (sweepindex==d-1 && ltr) || (sweepindex==1 && ~ltr)          % check whether half a sweep passed
        itr=itr+1;
        yhat=sim_volterraTN(u,TN);
        yhat=reshape(yhat',[N*l,1]);
        e(itr)=sqrt(norm(y-yhat)^2/N);    % RMS metric for residual
    end     
    updatesweep;    
end  

    function updateTT
        A=zeros(N*l,(p*M+1)^2*prod(r([sweepindex,sweepindex+2])));
        for i=M:N+M-1
            ui=zeros(p*M+1,1);
            ui(1)=1;
            for j=1:M
                ui(2+(j-1)*p:2+j*p-1)=uextended(i-j+1,:)';                
            end             
            vk1=eye(l);             % initialize row vector v_{k-1}
            vk2=1;                  % initialize column vector v_{k+2}
            for j=1:sweepindex-1
               vk1=vk1*reshape(reshape(permute(TN{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end
            for j=sweepindex+2:d
                vk2=vk2* reshape(reshape(permute(TN{j},[3 1 2]),[r(j+1)*r(j),p*M+1])*ui,[r(j+1),r(j)])';
            end        
            A((i-M)*l+1:(i-M+1)*l,:)=mkron(vk2',mkron(ui',2),vk1);
        end
        g=pinv(A)*y;                                        % truncated SVD solution (pseudoinverse)
        g=reshape(g,[r(sweepindex)*(p*M+1),(p*M+1)*r(sweepindex+2)]);
        [Ut,St,Vt]=svd(g);
        tol=eps(St(1))*max(size(g));
        st=diag(St);
        rankg=sum(st > tol);
        % we only update the TN rank such that update of all remaining
        % cores is guaranteed
        r2=r;
        r2(sweepindex+1)=rankg;
        rankg=checkRank(r2,N*l,(p*M+1)^2,sweepindex+1);
        r(sweepindex+1)=rankg;
        if ltr
            % left-to-right sweep, generate left orthogonal cores
            TN{sweepindex}=reshape(Ut(:,1:rankg),[r(sweepindex),(p*M+1),rankg]);
            TN{sweepindex+1}=reshape(St(1:rankg,1:rankg)*Vt(:,1:rankg)',[rankg,p*M+1,r(sweepindex+2)]);
        else
            % right-to-left sweep, generate right orthogonal cores
            TN{sweepindex}=reshape(Ut(:,1:rankg)*St(1:rankg,1:rankg),[r(sweepindex),(p*M+1),rankg]);
            TN{sweepindex+1}=reshape(Vt(:,1:rankg)',[rankg,p*M+1,r(sweepindex+2)]);
        end        
    end

    function newrank=checkRank(oldrank,Nsamples,nsquared,index)
        temp=zeros(1,length(oldrank)-2);
        for i=1:length(temp);
            temp(i)=oldrank(i)*nsquared*oldrank(i+2)-Nsamples;
        end
        if sum(temp < 0)==length(temp)
            newrank=oldrank(index);
        else
            [Y,I]=max(temp);
            temp=zeros(1,length(oldrank));
            temp(I)=oldrank(I);
            temp(I+2)=oldrank(I+2);
            temp(index)=0;
            newrank=floor(Nsamples/(nsquared*sum(temp)));
        end
    end


    function updatesweep
        if ltr
            if sweepindex < d-1
                sweepindex=sweepindex+1;
            else % sweepindex has reached end of the TN
                ltr=0; 
            end
        else
            if sweepindex > 1
                sweepindex=sweepindex-1;
            else % sweepindex has reached beginning of the TN
                ltr=1;
            end
        end
    end
end