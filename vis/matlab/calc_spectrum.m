function [range, spec_arr]=calc_spectrum(in,L,compensated)
%
% Plots magnetic field spectrum. 
% Tested for 3 dimensions. 
% Should work for any dimension
%
% Also will plot relevant wave-numbers 
% given in 'Simulations of the Small-Scale Turbulent Dynamo',
% Schekochihin et. al. 2004.
%

dim = ndims(in) - numel(find(size(in)==1));
ncell=max(size(in));

dk = 2 * pi / L;
idk=1.0/dk;
kmax = (ncell/2 -1)*dk;
kmin = - (ncell/2 )*dk;
[kx,ky,kz] = ndgrid(kmin:dk:kmax);
kx = ifftshift(kx);
ky = ifftshift(ky);
kz = ifftshift(kz);

ksqr = kx.^2 + ky.^2 + kz.^2;
k=sqrt(ksqr);
k= round(k*idk) + 1;

clearvars kx ky kz ksqr;

spec_arr=accumarray(k(:),in(:));
counts=accumarray(k(:),1);

clearvars k in_arr

range=0:dk:kmax; % anything larger than kmax is potentially underresolved

if(compensated) 
    range=0:dk:sqrt(dim)*kmax;
    l=length(range);
    
    surface = 2*pi^(dim/2.0)/(gamma(dim/2)*dk^(dim-1));
    comp=surface*range.^(dim-1)./(counts(1:l).');

    spec_arr=spec_arr(1:l).*comp.';
end
spec_arr = spec_arr*idk; %discrete to continuous

end
