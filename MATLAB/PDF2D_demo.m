% PDF2D_demo_class.m
% 
% Script used to estimate the joint cumulative density funtion (CDF) 
% and the joint probability density function (PDF) of a couple of random 
% variables by using a LUT interpolated by bidimensional cubic splines.
%
% This demo file implements the experimental tests shown in:
% M. Scarpiniti, R. Parisi and A. Uncini, "Flexible estimation of joint 
% probability and joint cumulative density functions", Electronics Letters, 
% Vol. 46, No. 15, p. 1084-1086, July 2010.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 19-Jan-2010$
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

disp('PDF2D_demo_class');

% Parameters --------------------------------------------------------------
aftype  = 2;        % 1: B-spline, 2: CR-spline
DeltaX  = 0.05;     % DeltaX
x_range = 2.1;      % Range limit
delta = 1e-5;       % Regularization factor
mu = 1e-7;          % Learning rate for ctrl points
ep = 10;            % Epochs of training

% Signal generation -------------------------------------------------------
M = 8;
L = 10000;
SNR = 20;
data = randi([0, M-1], L, 1);
txSig = pskmod(data, M, pi/M);  % M-PSK
% txSig = qammod(data, M);      % M-QAM
% x = awgn(txSig, SNR);
% x = randn(2, L)/6;
% x = rescale(x);
% x = rescale(x, -1, 1);
load 8PSK

% Nonlinearity definition and adaptation ----------------------------------
sp = LUT2D(x_range, DeltaX, aftype, mu, delta, ep);
sp.adapt(x);
sp.compute_CDF_PDF();
sp.make_plot();
% -------------------------------------------------------------------------
