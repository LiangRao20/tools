function [coeff,conf,dof] = dcregress(x, y, test_flag)

    if ~exist('test_flag', 'var'), test_flag = 0; end

    if test_flag
        test_dcregress;
        return;
    end

    if size(x,2) == 1, x = x'; end
    if size(y,2) == 1, y = y'; end

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

    dof(1) = calcdof(x);
    dof(2) = calcdof(y);
    dof = min(dof-2);

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
