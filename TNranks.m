function r=TNranks(TN)
% r=TNranks(TN)
% -------------
% Returns a vector containing the TN-ranks of a given MIMO Volterra tensor
% in the Tensor Network format TN.
%
% r         =   vector, r(i) contains TN-rank r_{i-1},
%
% TN        =   cell, TN{i} contains the ith TN core of the MIMO Volterra
%               tensor,
%
% Reference
% ---------
%
% 07-2016, Kim Batselier

d=length(TN);
r=ones(1,d+1);
for i=1:d
    r(i)=size(TN{i},1);
end

end