% Returns start and end inidices of gaps in input time series.
%    [gapstart, gapend] = FindGaps(in)
function [gapstart, gapend] = FindGaps(in)
    nans = isnan(in);
    edges = diff(nans);
    gapstart = find(edges == 1) + 1;
    gapend = find(edges == -1);

    if isnan(in(1)), gapstart = [1 gapstart]; end
    if isnan(in(end))
        gapend(end+1) = length(in);
    else
        % last gap is not at the end of time series.
        gapstart(end+1) = length(in);
        gapend(end+1) = length(in);
    end

    if isempty(gapstart) & isempty(gapend) ...
            & isequal(nans, zeros(size(nans)))
        % input series has no gaps
        gapstart = length(in) + 1;
        gapend = gapstart;
    end

    assert(length(gapstart) == length(gapend), ...
           ['FindGaps: gapstart and gapend are not same size.']);
end