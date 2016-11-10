function [YY_avg,freq]=spectrum_band_avg(yy,dt,M,winstr,plotflag);
%[spectrum,freq]=spectrum_band_avg(y,dt,M,winstr,plotflag);
%
%  Make one-sided, band-averaged spectral estimate for a real, 1-D data record
%
% Inputs:
%	y, vector to be transformed (must be real and uniformly spaced)
%	dt, time increment
%	M, number of bands to average
%	winstr, taper window; can be 'hann', 'blackman', 'parzenwin', 
%			'triang', 'gausswin', 'no taper', 
%			or 'end taper' (cosine taper on first and last 10% of segment)
%	plotflag, plots spectrum if set to 1
%
% Outputs:
%	spectrum, estimated spectrum
%	freq, frequency vector (in units of cycles/unit time)
%
%
% This code was written for MIT 12.805 class example.  Calls centeredFFT.m,
%  band_avg.m, confid.m
%
% Tom Farrar, 2016, jfarrar@whoi.edu


yy=yy(:);
N=length(yy);

if strcmp(winstr,'end taper')%Window with tapered cosine over first/last 10% of record
    N10=fix(N/5);
    win=ones(N,1);
    win(1:N10,1)=1-cos([1:N10]'*pi/(2*N10));
    win(N:-1:N-N10+1,1)=win(1:N10,1);
elseif strcmp(winstr,'hann')%use Hanning window and apply variance correction
    win=hann(N);win=win./sqrt(sum(win.^2.*1./N));
elseif strcmp(winstr,'blackman')%use Blackman window and apply variance correction
    win=blackman(N);win=win./sqrt(sum(win.^2.*1./N));
elseif strcmp(winstr,'parzenwin')%use Parzen window and apply variance correction
    win=parzenwin(N);win=win./sqrt(sum(win.^2.*1./N));
elseif strcmp(winstr,'triang')%use triangular window and apply variance correction
    win=triang(N);win=win./sqrt(sum(win.^2.*1./N));
elseif strcmp(winstr,'gausswin')%use gaussian window with alpha=2.5, and apply variance correction
    win=gausswin(N,3.5);win=win./sqrt(sum(win.^2.*1./N));
elseif strcmp(winstr,'no taper') %do nothing (ie, multiply by ones)
    win=ones(N,1);
    %display('Using rectangular window...')
end

T=N.*dt;
yy=yy-mean(yy);
yy=win.*yy;
[Y,freq_i]=centeredFFT(yy,dt);
ff=find(freq_i>0);
Y=Y(ff);
freq_i=freq_i(ff);
YY_raw=2.*T./N.^2.*Y.*conj(Y);

[YY_avg]=band_avg(YY_raw,M);
freq=band_avg(freq_i,M);

if plotflag==1
    [low,up]=confid(0.05,2*M);

    figure
    h1=loglog(freq,YY_avg)
    hold on%JTF
           %loglog(freq,fitspec,'r')
    spot=freq(4)*2;wid=spot*0.05;cstr=get(h1,'color');
    lowzz=var(yy)*low./100;
    upzz=var(yy)*up./100;
    if 1==1%plot errorbar for tapered case with EDoF
        loglog(spot,var(yy)./100,'o','color',cstr);
        h3=loglog([spot+10^-12 spot-10^-12], [lowzz,upzz],'color',cstr);
        loglog([spot+wid spot-wid], [upzz,upzz],'color',cstr);
        loglog([spot+wid spot-wid], [lowzz,lowzz],'color',cstr);
        hh=text(1.1*spot,var(yy)./100,'95%');set(hh,'horizontalalignment','left')
    end%if plot errorbar
    axis tight;ax=axis;axis([ax(1:3) ax(4).*1.05]);ax=axis;
    xlabel('Frequency (cyc/day)')
    ylabel('Spectral density (m^2/s^2/cpd)')
    %legend([h1 h2],'No taper',['With taper'])
    text(spot,10^-5,['M=' num2str(M)])
end









