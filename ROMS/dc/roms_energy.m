% function [EKE,MKE,PE] = roms_energy(fname,tindices)

% CHANGELOG
% First actual version                                          23 Feb 2012

% Todo list
% 1) add netcdf output
% 2) Is there a bug in PE calculation?
% 3) Try lower values of slab to see which works best
% 4) Add w contribution to KE

function [EKE,MKE,PE] = roms_energy(fname,tindices)

% input
% fname = 'his';
% fname = find_file(fname);
% tindices = [1 Inf];

% parameters
vinfo = ncinfo(fname,'u');
dim   = length(vinfo.Size); 
slab  = 40;

warning off
grid = roms_get_grid(fname,fname,0,1);
warning on

% parse input
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

%% read data

% caps indicates domain integrated values
EKE = nan(nt,1);
MKE = EKE;
PE  = EKE;

R0  = ncread(fname,'R0');
time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;

zrho = permute(grid.z_r,[3 2 1]);
try
    cpb = progressbar();
catch ME
    cpb = [];
end

for i=0:iend-1
    % FROM mod_movie.m - propagate changes back
    % start and count arrays for ncread : corrected to account for stride
    
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    
    if isempty(cpb), fprintf('\nReading Data...\n'); end
    u   = ncread(fname,'u',read_start,read_count,stride); pbar(cpb,i+1,1,iend,4);
    v   = ncread(fname,'v',read_start,read_count,stride); pbar(cpb,i+1,2,iend,4);
    %w   = ncread(fname,'w',read_start,read_count,stride); pbar(cpb,i+1,3,iend,4);
    rho = R0 + ncread(fname,'rho',read_start,read_count,stride); pbar(cpb,i+1,4,iend,4);
	if isempty(cpb), fprintf('\n Done reading data... \n'); end
    
    % mean fields
    um = mean(u,2);
    vm = mean(v,2);
    %wm = mean(w,2);
    rm = mean(rho,2);

    % eddy fields
    up = bsxfun(@minus,u,um);
    vp = bsxfun(@minus,v,vm);
    %wp = bsxfun(@minus,w,wm);
    rp = bsxfun(@minus,rho,rm);
    
    % average so that everything lands up on interior-rho points
    up = (up(1:end-1,2:end-1,:,:) + up(2:end,2:end-1,:,:))/2;
    um = (um(1:end-1, :     ,:,:) + um(2:end, :     ,:,:))/2;
    vp = (vp(2:end-1,1:end-1,:,:) + vp(2:end-1,2:end,:,:))/2;
    %vm = vm(2:end-1,:,:,:);
    
    eke = 0.5*rho(2:end-1,2:end-1,:,:).*(up.^2 + vp.^2); % SLOW?!
    mke = 0.5*bsxfun(@times,rho(2:end-1,2:end-1,:,:),(um.^2 + vm(2:end-1,:,:,:).^2));
    oke = rho(2:end-1,2:end-1,:,:).*(bsxfun(@times,up,um)+ bsxfun(@times,vp,vm(2:end-1,:,:,:)));
    pe  = 9.81*bsxfun(@times,rho,zrho);
    
    s = size(u);
    
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end);

    EKE(tstart:tend-1) = domain_integrate(eke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
    MKE(tstart:tend-1) = domain_integrate(mke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
    OKE(tstart:tend-1) = domain_integrate(oke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
    PE(tstart:tend-1)  = domain_integrate(pe(2:end-1,2:end-1,:,:),grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));  
end

if ~isempty(cpb)
    cpb.stop();
end

%% plot
figure;
plot(time,PE/nanmax(PE(:)),'b'); hold on;
plot(time,EKE/nanmax(EKE(:)),'r');
plot(time,MKE/nanmax(MKE(:)),'k');
plot(time,OKE/nanmax(OKE(:)),'m');
ylabel('Normalized Energy');
xlabel('Time (days)');
legend('PE','EKE','MKE','OKE');

figure;
hold on;
plot(time,EKE,'r');
plot(time,MKE,'k');
plot(time,OKE,'m');
ylabel('Energy');
xlabel('Time (days)');
legend('EKE','MKE','OKE');

figure;
plot(time,PE);
ylabel('Energy');
xlabel('Time (days)');
legend('PE');

save energy.mat time PE EKE MKE OKE

%% local functions

function [] = pbar(cpb,i,j,imax,jmax)
    if ~isempty(cpb)
        txt = sprintf(' Progress: i=%d, j=%d',i,j);
        progressbarupdate(cpb,(jmax*(i-1)+j)/(imax*jmax)*100,txt);
    end