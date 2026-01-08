# -*- coding: utf-8 -*-
"""
Script used to estimate the comulative density funtion (CDF) and the 
probability density function (PDF) of a random vector by using a LUT 
interpolated by cubic splines.

This demo file implements the experimental tests shown in:
- M. Scarpiniti, R. Parisi and A. Uncini, "Flexible estimation of probability 
and cumulative density functions", Electronics Letters, Vol. 45, N. 21, 2009.

Info: michele.scarpiniti@uniroma1.it
-------------------------------------------------------------------------------
$Revision: 1.0$  $Date: 28-Feb-2008$
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
from scipy import stats
from LUT import LUT  # Ensure the previous Python code is saved as LUT.py


def main():
    print('PDF_demo')

    # Parameters --------------------------------------------------------------
    aftype = 1        # 1: B-spline, 2: CR-spline
    delta_x = 0.05    # DeltaX
    x_range = 1.1     # Range limit
    delta = 1e-4      # Regularization factor
    mu = 1e-4         # Learning rate for ctrl points
    ep = 10           # Epochs of training

    # Signal generation -------------------------------------------------------
    L = int(1e5)
    
    # MATLAB: random('Exponential', 2, 1, L) 
    # Note: Scipy's 'scale' is the mean, same as MATLAB's second argument.
    x = stats.expon.rvs(scale=2, size=L)

    # Other signal options (commented out as in original):
    # x = stats.norm.rvs(loc=0, scale=5, size=L)
    # x = stats.uniform.rvs(loc=0, scale=1, size=L)

    
    # Rescale to [0, 1] as per the MATLAB rescale(x) default
    x_min, x_max = x.min(), x.max()
    x = (x - x_min) / (x_max - x_min)


    # Nonlinearity definition and adaptation ----------------------------------
    # Initialize the LUT object
    sp = LUT(x_range, delta_x, aftype, mu, delta, ep)
    
    print("Starting adaptation...")
    sp.adapt(x)              # Perform adaptation
    
    print("Computing CDF and PDF...")
    sp.compute_cdf_pdf()     # Estimate functions
    
    print("Generating plots...")
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
