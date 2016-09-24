% estimate degrees of freedom using integral timescale
%        [dof, IntegralTimeScale] = calcdof(in)
function [dof, IntegralTimeScale] = calcdof(in)

    [c,lags] = xcorr(in - mean(in), 'coef');
    [cmax,imax] = max(c);
    % i0 = imax + find(c(imax:end) < 0, 1, 'first') - 1; % first zero crossing

    % calculate a bunch of timescales and take maximum
    % From Talley et al., Descriptive Physical Oceanography.
    IntegralTimeScale = nan([1 length(c) - (imax + 3) + 1]);
    for tt=imax+3:length(c)
        ff = tt - (imax+3) + 1;
        CorrArea = sum( (c(imax:tt-1) + c(imax+1:tt))/2 );
        IntegralTimeScale(ff) = CorrArea/cmax;
    end

    IntegralTimeScale = max(IntegralTimeScale);
    dof = floor(length(in)/IntegralTimeScale);
end
