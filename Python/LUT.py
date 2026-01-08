# -*- coding: utf-8 -*-
"""
CLASS used to estimate the comulative density funtion (CDF) and the 
probability density function (PDF) of a random  vector by using a LUT 
interpolated by cubic splines. 

The class is used for implementing the experimental tests shown in:
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
import matplotlib.pyplot as plt


class LUT:
    """
    Class used to estimate the Cumulative Density Function (CDF) 
    and the Probability Density Function (PDF) of a random 
    vector by using a LUT interpolated by cubic splines.
    """

    def __init__(self, range_val, delta_x, basis, mu, delta, epochs):
        self.range = range_val
        self.delta_x = delta_x
        self.basis = basis
        self.mu = mu
        self.delta = delta
        self.epochs = epochs

        self.nq = int(np.ceil(2 * self.range / self.delta_x) + 1)
        
        # Private-like attributes
        self._q = None
        self._m = None
        self._t = np.zeros(4)
        self._dt = np.zeros(4)
        self._ddt = np.zeros(4)
        self._u = 0.0
        self._i = 0  # Python uses 0-based indexing
        self._kk = 5000
        self._x = None
        self._cdf = None
        self._pdf = None

        self._init_q()
        self._init_sp_matrix()
        

    def _init_q(self):
        # Generates the control points
        self._q = np.arange(-self.range, self.range + self.delta_x, self.delta_x)
        

    def _init_sp_matrix(self):
        if self.basis == 1:  # B-spline
            self._m = (1/6) * np.array([
                [-1,  3, -3,  1],
                [ 3, -6,  3,  0],
                [-3,  0,  3,  0],
                [ 1,  4,  1,  0]
            ])
        elif self.basis == 2:  # CR-spline
            self._m = 0.5 * np.array([
                [-1,  3, -3,  1],
                [ 2, -5,  4, -1],
                [-1,  0,  1,  0],
                [ 0,  2,  0,  0]
            ])
        else:
            raise ValueError("Spline basis should be 1 or 2!")
            

    def process(self, x):
        su = x / self.delta_x + (self.nq - 1) / 2
        self._i = int(np.floor(su))
        self._u = su - self._i

        # Boundary checks (Adjusted for 0-based indexing)
        if self._i < 0:
            self._i = 0
        if self._i > (self.nq - 4):
            self._i = self.nq - 4

        self._t = np.array([self._u**3, self._u**2, self._u, 1])
        self._dt = np.array([3 * self._u**2, 2 * self._u, 1, 0])
        self._ddt = np.array([6 * self._u, 2, 0, 0])
        
        # Matrix multiplication: T * M * Q_slice
        return self._t @ self._m @ self._q[self._i : self._i + 4]
    

    def adapt(self, x):
        self._x = np.array(x)
        l_len = len(self._x)
        
        for e in range(1, self.epochs + 1):
            # Cosine annealing for learning rate
            mu_n = 0.5 * self.mu * (1 + np.cos((e - 1) * np.pi / self.epochs))
            
            for n in range(l_len):
                y = self.process(self._x[n])
                idx = slice(self._i, self._i + 4)
                
                # Gradient update
                grad = self._dt @ self._m
                denom = (grad @ self._q[idx]) + self.delta
                
                self._q[idx] = self._q[idx] + mu_n * grad / denom
                
                # Constraints
                self._q = np.sort(self._q)
                self._q[0] = 0.0
                self._q[-1] = 1.0
        return self
    

    def compute_cdf_pdf(self, kk=5000):
        self._kk = kk
        x_lim = round(self.range)
        
        self._cdf = np.zeros(self._kk)
        self._pdf = np.zeros(self._kk)
        xa = np.zeros(self._kk)
        
        dx = 2 * x_lim / self._kk
        xx = -x_lim
        
        for k in range(self._kk):
            self._cdf[k] = self.process(xx)
            idx = slice(self._i, self._i + 4)
            self._pdf[k] = (self._dt @ self._m @ self._q[idx]) / self.delta_x
            xa[k] = xx
            xx += dx
            
        # Normalize PDF
        area = dx * np.sum(self._pdf)
        if area > 0:
            self._pdf = self._pdf / area
        return self

    
    def make_plot(self):
        if self._x is None or self._cdf is None:
            print("Data not processed. Run adapt() and compute_cdf_pdf() first.")
            return

        x_lim = round(self.range)
        xa = np.linspace(-x_lim, x_lim, self._kk)

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10))
        
        # Input Data
        ax1.plot(self._x, color='black')
        ax1.set_title('Input data')
        ax1.grid(True)
        ax1.set_ylabel('Value')
        ax1.set_xlabel('n')

        # CDF
        ax2.plot(xa, self._cdf, color='black', linewidth=2)
        ax2.set_title('Estimated CDF')
        ax2.set_xlim([-x_lim - 0.2, x_lim + 0.2])
        ax2.set_ylim([-0.2, 1.2])
        ax2.grid(True)
        ax2.set_ylabel('CDF')
        ax2.set_xlabel('x')

        # PDF
        ax3.plot(xa, self._pdf, color='black', linewidth=2)
        ax3.set_title('Estimated PDF')
        ax3.set_xlim([-x_lim - 0.2, x_lim + 0.2])
        ax3.grid(True)
        ax3.set_ylabel('PDF')
        ax3.set_xlabel('x')

        plt.tight_layout()
        plt.show()

   
    # Getters
    def get_q(self): return self._q
    def get_cdf(self): return self._cdf
    def get_pdf(self): return self._pdf


"""
-------------------------------------------------------------------------------
© 2025
M. Scarpiniti, 
Department of Information Engineering, Electronics and Telecommunications
(DIET) -- 'Sapienza' University of Rome
-------------------------------------------------------------------------------
"""
