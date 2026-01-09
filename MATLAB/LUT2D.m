% LUT2D.m
% 
% CLASS used to estimate the joint cumulative density funtion (CDF) 
% and the joint probability density function (PDF) of a couple of random 
% variables by using a LUT interpolated by bidimensional cubic splines.
%
% The class is used for implementing the experimental tests shown in:
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

classdef LUT2D < handle

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

        function obj = LUT2D(range, DeltaX, basis, mu, delta, epochs)
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
            % obj.Q = ((-obj.range:obj.DeltaX:obj.range).' + ...
            %     (-obj.range:obj.DeltaX:obj.range));
            [X, Y] = meshgrid(-obj.range:obj.DeltaX:obj.range);
            A = 0.5*(ones(size(X)) + tanh(X));
            B = 0.5*(ones(size(Y)) + tanh(Y));
            obj.Q = sqrt((A.^2 + B.^2)/2);
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
 
            obj.i(obj.i<1) = 1;                    % The index must start from 1

            obj.i(obj.i>(obj.Nq-3)) = obj.Nq - 3;  % The index cannot exceed np - 3

            ii1 = obj.i(1) : obj.i(1)+3;
            ii2 = obj.i(2) : obj.i(2)+3;
            obj.T   = [obj.u.^3, obj.u.^2, obj.u, ones(size(obj.u))];   % u^T
            obj.dT  = [3*obj.u.^2, 2*obj.u, ones(size(obj.u)), zeros(size(obj.u))];     % dot(u)^T
            obj.ddT = [6*obj.u, 2*ones(size(obj.u)), zeros(size(obj.u)), zeros(size(obj.u))];             % dot(dot(u))^T
            y = obj.T(2,:) * obj.M * (obj.T(1,:) * obj.M * obj.Q(ii1,ii2)).';   % Eq. (5): T2 M (T1 M Q)^T
        end

        function obj = adapt(obj, x)
            if size(x,2) < size(x,1)
                x = x.';
            end
            if ~isreal(x)
                re = real(x);
                im = imag(x);
                x = [re; im];
            end

            obj.x = x;
            L = length(x);
            for e = 1:obj.epochs
                mu_n = 0.5*obj.mu * (1 + cos((e-1)*pi/obj.epochs));  % Cosine annealing
                for n = 1:L
                    y = obj.process(x(:,n));
                    ii1 = obj.i(1) : obj.i(1)+3;
                    ii2 = obj.i(2) : obj.i(2)+3;
                    % Eq. (6)
                    num = abs(diag((obj.T(2,:) * obj.M) .* (obj.dT(1,:) * obj.M) - ...
                        (obj.dT(2,:) * obj.M) .* (obj.T(1,:) * obj.M)));
                    den = abs(obj.T(2,:) * obj.M * (obj.dT(1,:) * obj.M * obj.Q(ii1,ii2)).' - ...
                        obj.dT(2,:) * obj.M * (obj.T(1,:) * obj.M * obj.Q(ii1,ii2)).');
                    obj.Q(ii1,ii2) = obj.Q(ii1,ii2) + mu_n * num / (den + obj.delta);  % Eq. (7)
                
                    % Sort the spline
                    obj.Q = sort(sort(obj.Q, 2));
                    obj.Q(1:2, 1:2) = zeros(2);
                    obj.Q(end-1:end, end-1:end) = ones(2);
                end
            end
        end

        function obj = compute_CDF_PDF(obj, KK)
            % Compute the estimated CDF and PDF
            xLIM = round(obj.range);
            if nargin < 2
                obj.KK = 50;
            else
                obj.KK = KK;
            end
            obj.cdf = zeros(obj.KK, obj.KK);
            obj.pdf = zeros(obj.KK, obj.KK);
            dx = 2*xLIM/obj.KK;
            xx1 = -xLIM;
            for k1 = 1:obj.KK
                xx2 = -xLIM;
                for k2 = 1:obj.KK
                    obj.cdf(k1,k2) = obj.process([xx1, xx2].');   % CDF 
                    ii1 = obj.i(1) : obj.i(1)+3;
                    ii2 = obj.i(2) : obj.i(2)+3;
                    obj.pdf(k1,k2) = obj.dT(2,:) *obj.M * (obj.dT(1,:) * obj.M * obj.Q(ii1,ii2)).' / obj.DeltaX^2;  % PDF
                    % obj.pdf(k1,k2) = 0.5*abs(obj.T(2,:) *obj.M * (obj.dT(1,:) * obj.M * obj.Q(ii1,ii2)).' - ...
                    %     obj.dT(2,:) *obj.M * (obj.T(1,:) * obj.M * obj.Q(ii1,ii2)).')/ obj.DeltaX;  % PDF
                    xx2 = xx2 + dx;
                end
                xx1 = xx1 + dx;
            end

            % Remove border effects
            obj.pdf(1:2,:) = 0;
            obj.pdf(end-1:end,:) = 0;
            obj.pdf(:,1:2) = 0;
            obj.pdf(:,end-1:end) = 0;
            
            % Normalize the PDF
            volume = trapz(dx, trapz(dx, abs(obj.pdf), 2), 1);
            obj.pdf = abs(obj.pdf)/volume;
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

        function UI = get_ui(obj)
        % This function returns the local spline index i and abscissa u.
            UI = [obj.u, obj.i];
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
            xa = -xLIM:2*xLIM/obj.KK:xLIM;
            xa(end) = [];
            figure('PaperSize',[20.98 29.68]);
            box('on');

            % Plot the input data
            subplot(1,3,1);
            scatter(obj.x(1,:), obj.x(2,:), '.k');
            grid on;
            xlim([-max(abs(obj.x(1,:)))-0.2, max(abs(obj.x(1,:)))+0.2]);
            ylim([-max(abs(obj.x(2,:)))-0.2, max(abs(obj.x(2,:)))+0.2]);
            title('Input data');
            xlabel('Input {\itx}');
            ylabel('Input {\ity}');
            set(gca, 'FontSize', 10, 'FontWeight', 'demi');
            
            % Plot the estimated CDF
            subplot(1,3,2);
            mesh(xa, xa, obj.cdf, 'EdgeColor', 'k');
            grid on;
            xlim([-xLIM-0.2 xLIM+0.2]);
            ylim([-xLIM-0.2 xLIM+0.2]);
            zlim([-0.2 1.2]);
            title('Estimated CDF'); 
            xlabel('Input {\itx}');
            ylabel('Input {\ity}');
            zlabel('CDF');
            set(gca, 'FontSize', 10, 'FontWeight', 'demi');
            
            % Plot the estimated PDF
            subplot(1,3,3);
            mesh(xa, xa, abs(obj.pdf), 'EdgeColor', 'k');
            grid on;
            xlim([-xLIM-0.2 xLIM+0.2]);
            ylim([-xLIM-0.2 xLIM+0.2]);
            zlim([-0.2 max(obj.pdf(:))+0.2]);     
            title('Estimated PDF');
            xlabel('Input {\itx}');
            ylabel('Input {\ity}');
            zlabel('PDF');
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
