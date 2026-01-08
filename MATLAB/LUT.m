% LUT.m
% 
% CLASS used to estimate the cumulative density funtion (CDF) 
% and the probability density function (PDF) of a random 
% vector by using a LUT interpolated by cubic splines.
%
% The class is used for implementing the experimental tests shown in:
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

classdef LUT < handle

    properties
        range;
        DeltaX;
        basis;  % 1: B-spline, 2: CR-spline
        mu;
        delta;
        epochs;
    end

    properties (Access=private)
        Nq;
        Q;
        M;
        T = zeros(4, 1);
        dT = zeros(4, 1);
        ddT = zeros(4, 1);
        u = 0;
        i = 1;
        KK;
        x;
        cdf;
        pdf;
    end

    methods

        function obj = LUT(range, DeltaX, basis, mu, delta, epochs)
            obj.range = range;
            obj.DeltaX = DeltaX;
            obj.basis = basis;
            obj.mu = mu;
            obj.delta = delta;
            obj.epochs = epochs;

            obj.Nq = ceil(2*obj.range/obj.DeltaX) + 1;
            obj.init_Q();
            obj.init_sp_matrix();
        end

        function init_Q(obj)
            obj.Q = (-obj.range:obj.DeltaX:obj.range).';
        end

        function init_sp_matrix(obj)
            if obj.basis == 1
                obj.M = (1/6)*[-1  3 -3  1; ... 
                    3 -6  3  0; ... 
                    -3  0  3  0; ... 
                    1  4  1  0];
            elseif obj.basis == 2
                obj.M = 0.5*[-1  3 -3  1; ...
                    2 -5  4 -1; ...  
                    -1  0  1  0; ... 
                    0  2  0  0];
            else
                obj.M = zeros(4);
                error('Spline basis should be 1 or 2!');
            end
        end

        function y = process(obj, x)
            Su = x/obj.DeltaX + (obj.Nq-1)/2;
            obj.i = floor(Su);   % Local index
            obj.u = Su - obj.i;  % Local abscissa
 
            if obj.i<1                  % The index must start from 1
                obj.i = 1;
            end
            if obj.i>(obj.Nq-3)             % The index cannot exceed np - 3
                obj.i = obj.Nq - 3;
            end

            obj.T   = [obj.u^3 obj.u^2 obj.u 1];   % u^T
            obj.dT  = [3*obj.u^2 2*obj.u 1 0];     % dot(u)^T
            obj.ddT = [6*obj.u 2 0 0];             % dot(dot(u))^T
            y = obj.T * obj.M * obj.Q(obj.i : obj.i+3);   % Eq. (5): T M Q_i
        end

        function obj = adapt(obj, x)
            obj.x = x;
            L = length(x);
            for e = 1:obj.epochs
                mu_n = 0.5*obj.mu * (1 + cos((e-1)*pi/obj.epochs));  % Cosine annealing
                for n = 1:L
                    y = obj.process(x(n));
                    ii = obj.i : obj.i+3;
                    % obj.Q(ii) = obj.Q(ii) + obj.mu * (obj.dT*obj.M).' / ...
                    %     (obj.dT*obj.M*obj.Q(ii) + obj.delta);  % Eqs. (6)-(7)
                    obj.Q(ii) = obj.Q(ii) + mu_n * (obj.dT*obj.M).' / ...
                        (obj.dT*obj.M*obj.Q(ii) + obj.delta);  % Eqs. (6)-(7)
                
                    % Sort the spline
                    obj.Q = sort(obj.Q);
                    obj.Q(1) = 0;
                    obj.Q(end) = 1;
                end
            end
        end

        function obj = compute_CDF_PDF(obj, KK)
            % Compute the estimated CDF and PDF
            xLIM = round(obj.range);
            if nargin < 2
                obj.KK = 5000;
            else
                obj.KK = KK;
            end
            xa  = zeros(1, obj.KK);
            obj.cdf = zeros(1, obj.KK);
            obj.pdf = zeros(1, obj.KK);
            dx = 2*xLIM/obj.KK;
            xx = -xLIM;   
            for k = 1:obj.KK
                obj.cdf(k) = obj.process(xx);   % CDF 
                ii = obj.i : obj.i+3;
                obj.pdf(k) = obj.dT * obj.M * obj.Q(ii) / obj.DeltaX;  % PDF
                xa(k) = xx;
                xx = xx + dx;
            end
            
            % Normalize the PDF
            area = dx*sum(obj.pdf);
            obj.pdf = obj.pdf./area;
        end


        function q = get_Q(obj)
        % This function returns the spline control points.
            q = obj.Q;
        end

        function t = get_T(obj)
        % This function returns the T vector.
            t = obj.T;
        end

        function dt = get_dT(obj)
        % This function returns the dotT vector.
            dt = obj.dT;
        end

        function ddt = get_ddT(obj)
        % This function returns the dotdotT vector.
            ddt = obj.ddT;
        end

        function [u, i] = get_ui(obj)
        % This function returns the local spline index i and abscissa u.
            u = obj.u;
            i = obj.i;
        end

        function m = get_M(obj)
        % This function returns the spline basis matrix.
            m = obj.M;
        end

        function s = get_x(obj)
        % This function returns the input data.
            s = obj.x;
        end

        function f = get_cdf(obj)
        % This function returns the estimated CDF.
            f = obj.cdf;
        end

        function p = get_pdf(obj)
        % This function returns the estimated PDF.
            p = obj.pdf;
        end        

        function make_plot(obj)
            % Plot a single figure with data, CDF, and PDF
            xLIM = round(obj.range);
            xa = -xLIM:2*xLIM/5000:xLIM;
            xa(end) = [];
            figure('PaperSize',[20.98 29.68]);
            box('on');
            
            % Plot the input data
            subplot(3,1,1);
            plot(obj.x, 'k');
            grid on;
            ylim([min(obj.x), max(obj.x)+0.2]);
            title('Input data');
            xlabel('Index [{\itn}]');
            ylabel('Value');
            set(gca, 'FontSize', 10, 'FontWeight', 'demi');
            
            % Plot the estimated CDF
            subplot(3,1,2);
            plot(xa, obj.cdf, 'k', 'LineWidth', 2);
            grid on;
            xlim([-xLIM-0.2 xLIM+0.2]);
            ylim([-0.2 1.2]);
            title('Estimated CDF'); 
            xlabel('Input {\itx}');
            ylabel('CDF');
            set(gca, 'FontSize', 10, 'FontWeight', 'demi');
            
            % Plot the estimated PDF
            subplot(3,1,3);
            plot(xa, obj.pdf, 'k', 'LineWidth', 2);
            grid on;
            xlim([-xLIM-0.2 xLIM+0.2]);
            ylim([-0.2 max(obj.pdf)+0.2]);     
            title('Estimated PDF');
            xlabel('Input {\itx}');
            ylabel('PDF');
            set(gca, 'FontSize', 10, 'FontWeight', 'demi');
        end

    end

end

% -------------------------------------------------------------------------
% © 2025
% M. Scarpiniti, 
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------
