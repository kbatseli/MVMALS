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

r=[1 TN.n(:,end)'];

end