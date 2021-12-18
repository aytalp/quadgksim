# quadgksim
Adaptive Gauss-Kronrod quandrature for simultaneous integration of similar integrands.

Upgraded by Aytac Alparslan using quadgk routine of Matlab. 2021
Based on "quadva" by Lawrence F. Shampine.
Ref: L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab", Journal of Computational and Applied Mathematics 211, 2008, pp.131-140.

In order to test the routine, download the files "quadgk_sim.m", "test_func.m" and "compare_time.m" and run the following command in Matlab (R2021a or older):
> compare_time(1e-14,1000);

This will integrate the functions given in "test_func.m" 1000 times with the relative error of 1e-14 along the triangular integration contour [0]-->[1+1i]-->[1-1i]-->[0] on the complex z-plane and compare the average time requirements and the integration results.
