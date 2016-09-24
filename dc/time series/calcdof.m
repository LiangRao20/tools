% estimate degrees of freedom using integral timescale
%        [dof, IntegralTimeScale] = calcdof(in)
function [dof, IntegralTimeScale] = calcdof(in)

    [c,lags] = xcorr(in - mean(in), 'coef');
    [cmax,imax] = max(c);
    % i0 = imax + find(c(imax:end) < 0, 1, 'first') - 1; % first zero crossing

    % calculate a bunch of timescales and take maximum
    % From Talley et al., Descriptive Physical Oceanography.
    IntegralTimeScale = max(cumtrapz(c(imax:length(c))))./cmax;

    dof = floor(length(in)/IntegralTimeScale);
end
