# -*- coding: utf-8 -*-
"""
Script used to estimate the joint comulative density funtion (CDF) and the 
joint probability density function (PDF) of a couple of random variables by 
using a LUT interpolated by bidimensional cubic splines.

This demo file implements the experimental tests shown in:
- M. Scarpiniti, R. Parisi and A. Uncini, "Flexible estimation of joint 
probability and joint cumulative density functions", Electronics Letters, 
Vol. 46, No. 15, p. 1084-1086, July 2010.

Info: michele.scarpiniti@uniroma1.it
-------------------------------------------------------------------------------
$Revision: 1.0$  $Date: 19-Jan-2010$
$Revision: 1.1$  $Date: 28-Dec-2025$
-------------------------------------------------------------------------------
License to use and modify this code is granted freely without warranty to all, 
as long as the original authors are referenced and attributed as such.
The original authors maintain the right to be solely associated with this work.
-------------------------------------------------------------------------------
© 2025
M. Scarpiniti, 
Department of Information Engineering, Electronics and Telecommunications
(DIET) -- 'Sapienza' University of Rome
-------------------------------------------------------------------------------
"""

import numpy as np
from LUT2D import LUT2D  # Assumes the previous LUT2D code is saved as LUT2D.py


def main():
    print('PDF2D_demo')

    # Parameters --------------------------------------------------------------
    aftype = 2        # 1: B-spline, 2: CR-spline
    delta_x = 0.05    # DeltaX
    x_range = 2.1     # Range limit
    delta = 1e-5      # Regularization factor
    mu = 1e-7         # Learning rate for ctrl points
    ep = 10           # Epochs of training

    # Signal generation -------------------------------------------------------
    # Replicating the 8-PSK logic from the MATLAB script
    M = 8
    L = 10000
    
    # Generate random data symbols
    data = np.random.randint(0, M, L)
    
    # M-PSK Modulation: exp(j * (2*pi*data/M + pi/M))
    tx_sig = np.exp(1j * (2 * np.pi * data / M + np.pi / M))
    
    # Add noise to simulate a realistic scenario (similar to awgn in MATLAB)
    snr_db = 30
    snr_linear = 10**(snr_db / 10)
    noise_volts = np.sqrt(1 / (2 * snr_linear))
    noise = noise_volts * (np.random.randn(L) + 1j * np.random.randn(L))
    
    x_complex = tx_sig + noise
    
    # Format x as a 2xL matrix (Real and Imaginary parts) for LUT2D
    x = np.vstack([x_complex.real, x_complex.imag])

    # Note: In the MATLAB script, a 'load 8PSK' was used. 
    # Here we generated it dynamically above to ensure the code is runnable.

    # Nonlinearity definition and adaptation ----------------------------------
    # Initialize the 2D LUT object
    sp = LUT2D(x_range, delta_x, aftype, mu, delta, ep)
    
    print("Starting adaptation on 2D data...")
    sp.adapt(x)              # Perform adaptation
    
    print("Computing joint CDF and PDF...")
    sp.compute_cdf_pdf()     # Estimate functions
    
    print("Generating 3D plots...")
    sp.make_plot()           # Visualize results
    # -------------------------------------------------------------------------



if __name__ == "__main__":
    main()


"""
-------------------------------------------------------------------------------
© 2025
M. Scarpiniti, 
Department of Information Engineering, Electronics and Telecommunications
(DIET) -- 'Sapienza' University of Rome
-------------------------------------------------------------------------------
"""
