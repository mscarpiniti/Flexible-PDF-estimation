# -*- coding: utf-8 -*-
"""
CLASS used to estimate the joint comulative density funtion (CDF) and the joint 
probability density function (PDF) of a couple of random variables by using a 
LUT interpolated by bidimensional cubic splines.

The class is used for implementing the experimental tests shown in:
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
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D


class LUT2D:
    """
    Class used to estimate the joint Cumulative Density Function (CDF) 
    and the joint Probability Density Function (PDF) of a couple of random 
    variables by using a LUT interpolated by bidimensional cubic splines.
    """

    def __init__(self, range_val, delta_x, basis, mu, delta, epochs):
        self.range = range_val
        self.delta_x = delta_x
        self.basis = basis
        self.mu = mu
        self.delta = delta
        self.epochs = epochs

        self.nq = int(np.ceil(2 * self.range / self.delta_x) + 1)
        
        # Internal state
        self._q = None
        self._m = None
        self._t = np.zeros((2, 4))
        self._dt = np.zeros((2, 4))
        self._ddt = np.zeros((2, 4))
        self._u = np.zeros(2)
        self._i = np.zeros(2, dtype=int)
        self._kk = 50
        self._x = None
        self._cdf = None
        self._pdf = None

        self._init_q()
        self._init_sp_matrix()
        

    def _init_q(self):
        # Generate meshgrid for initialization
        vec = np.arange(-self.range, self.range + self.delta_x, self.delta_x)
        x_mesh, y_mesh = np.meshgrid(vec, vec)
        
        # Initializing Q based on tanh surfaces
        a = 0.5 * (1 + np.tanh(x_mesh))
        b = 0.5 * (1 + np.tanh(y_mesh))
        self._q = np.sqrt((a**2 + b**2) / 2)
        

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
        # x should be a vector of size 2
        su = x / self.delta_x + (self.nq - 1) / 2
        self._i = np.floor(su).astype(int)
        self._u = su - self._i

        # Boundary checks (0-indexed)
        self._i = np.clip(self._i, 0, self.nq - 4)

        # Build T, dT, ddT for both dimensions
        # Row 0: x dimension, Row 1: y dimension
        self._t = np.column_stack([self._u**3, self._u**2, self._u, np.ones(2)])
        self._dt = np.column_stack([3*self._u**2, 2*self._u, np.ones(2), np.zeros(2)])
        self._ddt = np.column_stack([6*self._u, 2*np.ones(2), np.zeros(2), np.zeros(2)])

        ii1 = slice(self._i[0], self._i[0] + 4)
        ii2 = slice(self._i[1], self._i[1] + 4)
        
        # Eq. (5): T2 M (T1 M Q)^T
        q_slice = self._q[ii1, ii2]
        inner = self._t[0, :] @ self._m @ q_slice.T
        y = self._t[1, :] @ self._m @ inner.T
        return y
    

    def adapt(self, x):
        # Ensure x is in the correct shape (2, L)
        x = np.atleast_2d(x)
        if x.shape[0] != 2:
            x = x.T
        
        # Handle complex data
        if np.iscomplexobj(x):
            x = np.vstack([x.real, x.imag])

        self._x = x
        l_len = x.shape[1]
        
        for e in range(1, self.epochs + 1):
            mu_n = 0.5 * self.mu * (1 + np.cos((e - 1) * np.pi / self.epochs))
            
            for n in range(l_len):
                y = self.process(x[:, n])
                ii1 = slice(self._i[0], self._i[0] + 4)
                ii2 = slice(self._i[1], self._i[1] + 4)
                
                # Adaptation Logic
                term1 = (self._t[1, :] @ self._m).reshape(-1, 1) @ (self._dt[0, :] @ self._m).reshape(1, -1)
                term2 = (self._dt[1, :] @ self._m).reshape(-1, 1) @ (self._t[0, :] @ self._m).reshape(1, -1)
                num = np.abs(term1 - term2)
                
                q_local = self._q[ii1, ii2]
                val1 = self._t[1, :] @ self._m @ (self._dt[0, :] @ self._m @ q_local.T).T
                val2 = self._dt[1, :] @ self._m @ (self._t[0, :] @ self._m @ q_local.T).T
                den = np.abs(val1 - val2)

                self._q[ii1, ii2] += mu_n * num / (den + self.delta)
                
                # Sorting and constraints
                self._q = np.sort(np.sort(self._q, axis=1), axis=0)
                self._q[0:2, 0:2] = 0.0
                self._q[-2:, -2:] = 1.0
        return self
    

    def compute_cdf_pdf(self, kk=50):
        self._kk = kk
        x_lim = round(self.range)
        self._cdf = np.zeros((self._kk, self._kk))
        self._pdf = np.zeros((self._kk, self._kk))
        
        dx = 2 * x_lim / self._kk
        xx1 = -x_lim
        
        for k1 in range(self._kk):
            xx2 = -x_lim
            for k2 in range(self._kk):
                self._cdf[k1, k2] = self.process(np.array([xx1, xx2]))
                
                ii1 = slice(self._i[0], self._i[0] + 4)
                ii2 = slice(self._i[1], self._i[1] + 4)
                
                # PDF calculation based on partial derivatives
                q_local = self._q[ii1, ii2]
                inner_pdf = self._dt[0, :] @ self._m @ q_local.T
                self._pdf[k1, k2] = (self._dt[1, :] @ self._m @ inner_pdf.T) / (self.delta_x**2)
                xx2 += dx
            xx1 += dx

        # Remove border effects
        self._pdf[0:2, :] = 0; self._pdf[-2:, :] = 0
        self._pdf[:, 0:2] = 0; self._pdf[:, -2:] = 0
        
        # Normalize joint PDF
        volume = np.trapz(np.trapz(np.abs(self._pdf), dx=dx, axis=1), dx=dx)
        if volume > 0:
            self._pdf = np.abs(self._pdf) / volume
        return self
    

    def make_plot(self):
        x_lim = round(self.range)
        xa = np.linspace(-x_lim, x_lim, self._kk)
        ya = np.linspace(-x_lim, x_lim, self._kk)
        X, Y = np.meshgrid(xa, ya)

        fig = plt.figure(figsize=(18, 5))
        
        # 1. Scatter plot of input data
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.scatter(self._x[0, :], self._x[1, :], s=1, c='k', alpha=0.5)
        ax1.set_title("Input Data")
        ax1.set_xlabel(r"$x_1$")
        ax1.set_ylabel(r"$x_2$")
        ax1.grid(True)

        # 2. Joint CDF Surface
        ax2 = fig.add_subplot(1, 3, 2, projection='3d')
        ax2.plot_surface(X, Y, self._cdf.T, cmap='Greys', edgecolor='k', alpha=0.8)
        ax2.set_title("Estimated Joint CDF")
        ax2.set_zlim([-0.2, 1.2])
        ax2.xaxis.set_inverted(True)
        ax2.set_xlabel(r"$x_1$")
        ax2.set_ylabel(r"$x_2$")
        # ax2.view_init(45,45)

        # 3. Joint PDF Surface
        ax3 = fig.add_subplot(1, 3, 3, projection='3d')
        ax3.plot_surface(X, Y, self._pdf.T, cmap='Greys', edgecolor='k', alpha=0.8)
        ax3.set_title("Estimated Joint PDF")
        ax3.xaxis.set_inverted(True)
        ax3.set_xlabel(r"$x_1$")
        ax3.set_ylabel(r"$x_2$")
        # ax3.view_init(45,45)
        
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
