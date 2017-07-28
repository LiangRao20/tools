function  [phii,Var] = oa1(x,a,b,V,E,xm,phi,era)

% Do 1-D objective analysis with correct mean removal
% call as
%	[phii,Var] = oa1(x,a,b,V,E,xm,phi,era)
% where
%	x = grid to map to
%	a = Gaussian scale
%	b = cosine scale
%	V = variance of true phi
%	E = error variance
%	xm = observation points
%	phi = values measured at xm
%	phii = values on grid x
%	Var = error variance in fit
%	era = maximum allowable allowable error variance (percent)

% covariance function = 
%	V*exp(-(delx/a)^2)*cos(delx/b)

xm(isnan(phi)) = [];
phi = cut_nan(phi);

[nm,mm] = size(phi);

if nm > mm
   phi = phi';
   [nm,mm] = size(phi);
end

A = corc(xm,a,b,V,E);
Ai = inv(A);

C = calcm(xm,x,a,b,V);

xq = Ai*phi';
xmm = sum(xq);
ymm = sum(sum(Ai));
phm = xmm/ymm;

phii = C*(Ai*(phi-phm*ones(nm,mm))');
[ni,mi] = size(phii);
phii = phii + phm*ones(ni,mi);

Cxx = calcm(x,x,a,b,V);
xq = 1-sum(C*Ai,2);

em  = (xq*xq')/sum(sum(Ai));
Var = Cxx -C*Ai*C' + em;        %error variance

IVar = diag(Var);
i = find((IVar/V) > (era/100));

[n,m] = size(i);
phii(i) = NaN*ones(n,m);