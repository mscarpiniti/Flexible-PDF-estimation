% PDF_demo_class.m
% 
% Script used to estimate the cumulative density funtion (CDF) 
% and the probability density function (PDF) of a random 
% vector by using a LUT interpolated by cubic splines.
%
% This demo file implements the experimental tests shown in:
% M. Scarpiniti, R. Parisi and A. Uncini, "Flexible estimation of
% probability and cumulative density functions", Electronics Letters, Vol.
% 45, N. 21, 2009.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 28-Feb-2008$
% $Revision: 1.1$  $Date: 28-Dec-2025$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2025
% M. Scarpiniti, 
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------

disp('PDF_demo_class');

% Parameters --------------------------------------------------------------
aftype  = 1;        % 1: B-spline, 2: CR-spline
DeltaX  = 0.05;     % DeltaX
x_range = 1.1;      % Range limit
delta = 1e-4;       % Regularization factor
mu = 1e-4;          % Learning rate for ctrl points
ep = 10;            % Epochs of training

% Signal generation -------------------------------------------------------
L = 1e5;
x = random('Exponential', 2, 1, L);
% x = random('Normal', 0, 5, 1, L);
% x = random('Uniform', 0, 1, 1, L);
x = rescale(x);
% x = rescale(x, -1, 1);

% Nonlinearity definition and adaptation ----------------------------------
sp = LUT(x_range, DeltaX, aftype, mu, delta, ep);
sp.adapt(x);
sp.compute_CDF_PDF();
sp.make_plot();
% -------------------------------------------------------------------------
