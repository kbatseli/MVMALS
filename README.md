MIMO Volterra Modified Alternating Linear Scheme (MVMALS) for Matlab&copy;/Octave&copy;
--------------------------------------------------------------------------------------------------

This package contains Matlab/Octave code for the identification of MIMO Volterra systems using the (Modified) Alternating Linear Scheme (MALS). Tensor Networks are implemented as structures with two fields: *core* and *n*. The *core* field contains a cell of TN-cores. The *n* field contains a matrix of TN-ranks, where size(n,1) is the number of cores and size(n,2) indicates the order of the TN-core. In order to maintain consistency the values of n(1,1) and n(end,end) always need to be 1.

1. Functions
------------

* demo

Demonstrates the usage of the two algorithms ALS and MALS in the identification of SISO and MIMO Volterra systems.

* yhat=sim_volterraTN(u,TN)

Simulates a truncated MIMO Volterra series in the Tensor Network (TN) format for given inputs u(:,1),u(:,2),....

* [TN,e]=mvals(y,u,M,r,THRESHOLD)

MIMO Volterra Alternating Linear Scheme (MVALS) algorithm for solving the MIMO Volterra system identification problem in the Tensor Network format.

* [TN,e]=mvmals(y,u,M,d,THRESHOLD)

MIMO Volterra Modified Alternating Linear Scheme (MVMALS) algorithm for solving the MIMO Volterra system identification problem in the Tensor Network format. Upper bounds on the TN-ranks are limited by the total number of output samples.

* U=makeU(u,M,d)

Generates the N x (pM+1)^d matrix U such that y=U\*vec(V) for a p-input Volterra system, where vec(V) is the vectorization of the MISO Volterra tensor and N=size(u,1) is the number of measured samples. For the MIMO case one needs to compute y=U\*[vec(V1) vec(V2) ... vec(Vl)] for l separate Volterra tensors.


2. Reference
------------

"Tensor Network alternating linear scheme for MIMO Volterra system identification"

Authors: Kim Batselier, Zhongming Chen, Ngai Wong
