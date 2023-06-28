"# Radar-Clutter-Nonlinear-Shrinkage-Approach" 
This code repository is for the paper titled as "Radar Clutter Covariance Estimation: A Nonlinear Spectral Shrinkage Approach" <https://arxiv.org/abs/2302.02045>.
The Repository uses another repostory from the paper "D. L. Donoho , M. Gavish and I. M. Johnstone, "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL", <https://arxiv.org/abs/1311.0851>, to use the function optimal_shrinkage.m.
SDR contains the the function optimal_shrinkage.m.
Challenge Dataset can be accessed from <https://www.islinc.com/products/rfview>.
To get SNR vs Datapoints plot run rho_bound.m.
To get SNR vs Angle plot, in the rho_bounds.m plots fix N1=1024, and take the average over Doppler.
To get SNR vs Doppler plot, in the rho_bounds.m plots fix N1=1024, and take the average over Angle.
To get the Pd vs SNR plot run main.m.
To get MSE vs SNR plot run MSE.m.
Contact me: <sj474 "at" cornell.edu> for any issues.
