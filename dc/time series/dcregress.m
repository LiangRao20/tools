% Does Ordinary least squares. Calculates confidence bounds after
% estimating degrees of freedom.
%    [coeff,conf,dof] = dcregress(x, y, dof, test_flag, plot_flag)
function [coeff,conf,dof] = dcregress(x, y, dof, test_flag, plot_flag)

    if ~exist('dof', 'var'), dof = NaN; end
    if ~exist('test_flag', 'var'), test_flag = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test_flag
        test_dcregress;
        return;
    end

    if size(x,2) == 1, x = x'; end
    if size(y,2) == 1, y = y'; end

    xnan = x; ynan = y;
    x = cut_nan(x); y = cut_nan(y);

    E = [ones(size(x')) x'];

    %%%%%%%%%%% See Wunsch(1996) pg. 116
    % P matrix
    coeff = E\y';
    intercept = coeff(1);
    slope = coeff(2);
    true = y; est = intercept + slope .* x;
    res = true-est;
    % (E' * E) ^-1
    %ETEI = inv(E'*E);
    % from http://blogs.mathworks.com/loren/2007/05/16/purpose-of-inv/
    [Q,R] = qr(E,0);
    S = inv(R);
    ETEI = S*S';
    % assuming noise vector (res) is white
    P = ETEI * E' * var(res) * E * ETEI;
    err = sqrt(diag(P));  % standard error

    if isnan(dof)
        % calculate degrees of freedom taking into account gaps
        dof(1) = calcdof(xnan);
        dof(2) = calcdof(ynan);
        dof = min(dof-2);
    end

    [lo,up] = conft(0.05, dof);
    conf(1,1) = up * err(1);
    conf(2,1) = up * err(2);

    %%%%%%%%%%% use MATLAB regress
    % [b, bint, r, rint, stats] = regress(y', E);

    % slope = b(2);
    % err = abs(bint(2) - b(2));

    % plot fit
    % figure; hold all;
    % plot(tvec/86400, true, '*');
    % plot(tvec/86400, est); plot(tvec/86400, res); liney(0);

    %[c,lags] = xcorr(fluxvec - mean(fluxvec(:)), 'coef');
    %plot(lags, c); linex(0); liney(0);

    if plot_flag
        figure; hold on;
        plot(x,y, '.', 'MarkerSize', 16, 'Color', [1 1 1]*0.4);
        limx = xlim;
        xvec = linspace(limx(1), limx(2), 1000);
        yvec = slope * xvec + intercept;
        yvechi = (slope + conf(2))*xvec + intercept + conf(1);
        yveclo = (slope - conf(2))*xvec + intercept - conf(1);
        plot(xvec, yvec, 'k');
        plot(xvec, yvechi, 'k--');
        plot(xvec, yveclo, 'k--');
    end
end

function [] = test_dcregress()

    b = 0.1
    a = 0.5
    x = linspace(0,1,8000);
    y = a*x + b + a/4 * rand(size(x));

    [coef,conf] = dcregress(x,y);
    [coef-conf coef+conf]
    % [b, bint, r, rint, stats] = regress(y', [ones(size(x')) x']);
end
