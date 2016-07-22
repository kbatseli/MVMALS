MIMO Volterra Modified Alternating Linear Scheme (MVMALS) for Matlab&copy;/Octave&copy;
--------------------------------------------------------------------------------------------------

This package contains Matlab/Octave code for the identification of MIMO Volterra systems using the (Modified) Alternating Linear Scheme (MALS). 

1. Functions
------------

* yhat=sim_volterraTT(u,TT)

Simulates a truncated MISO Volterra series in the Tensor Train (TT) format for given inputs u(:,1),u(:,2),....

* [TT,e]=mvals(y,u,M,r,THRESHOLD)

MIMO Volterra Alternating Linear Scheme (MVALS) algorithm for solving the MIMO Volterra system identification problem in the Tensor Train format. This function handles the MISO case, use it for each output separately.

* [TT,e]=mvmals(y,u,M,r,THRESHOLD)

MIMO Volterra Modified Alternating Linear Scheme (MVMALS) algorithm for solving the MIMO Volterra system identification problem in the Tensor Train format. This function handles the MISO case, use it for each output separately.


2. Reference
------------

"Tensor Train alternating linear scheme for MIMO Volterra system identification"

Authors: Kim Batselier, Zhongming Chen, Ngai Wong
