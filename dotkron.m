function y=dotkron(varargin)
% y=dotkron(A,B) or y=dotkron(A,B,C)
% ----------------------------------
% Computes the row-wise right-Kronecker product of matrices A,B,C. Note
% that this computes the right-Kronecker product such that the order of the
% indices is maintained.
%
% y,A,B,C   = matrix.
%
% Reference
% ---------
%
% 11/2016, Zhongming Chen

if nargin==2
    L=varargin{1};     R=varargin{2};
    [r1,c1]=size(L);   [r2,c2]=size(R);
    if r1 ~= r2
        error('Matrices should have equal rows!');
    else
        y=repmat(L,1,c2).*kron(R, ones(1, c1));  % faster than using FOR
    end
elseif nargin==3
    L=varargin{1};     M=varargin{2};     R=varargin{3};
    r1=size(L,1);      r2=size(M,1);      r3=size(R,1);
    if r1 ~= r2 || r2 ~=r3  ||  r1~=r3
        error('Matrices should have equal rows!');
    else
        y=dotkron(L, dotkron(M, R));
    end
else
    error('Please input 2 or 3 matrices!');
end
    


end