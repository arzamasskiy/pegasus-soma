function B_spectrum(U,basename_in)
%
% Plots magnetic field spectrum. 
% Tested for 3 dimensions. 
% Should work for any dimension
%
% Also will plot relevant wave-numbers 
% given in 'Simulations of the Small-Scale Turbulent Dynamo',
% Schekochihin et. al. 2004.
%

basename='/tigress/dstonge/pegasus/huge_run';
compensated=true; % try to compensate for grid cut-offs
if(nargin ==2)
   basename = basename_in;
end    

if(nargin ==0)
  % Read in velocity field   VTK
  fname = '/tigress/dstonge/pegasus/prod_run2/combined/combined.0150.fld.vtk'; % filename
  %fname = '/tigress/dstonge/pegasus/run7_alt2/combined/dynamo.0050.fld.vtk'; % filename

  % open file and initialize grid
  mom_ary = [3 3];
  [Grid,~] = init_grid(fname,mom_ary);

  U.L = Grid.x1max - Grid.x1min;

  for varid=1:length(mom_ary)
    [U.time,var,~,~] = readvtk(Grid,fname,varid);
    if (varid==1)        % read cell-centered magnetic field
      U.bx = squeeze(var(1,:,:,:));
      U.by = squeeze(var(2,:,:,:));
      U.bz = squeeze(var(3,:,:,:));
    elseif (varid==2)    % read cell-centered electric field
      U.ex = squeeze(var(1,:,:,:));
      U.ey = squeeze(var(2,:,:,:));
      U.ez = squeeze(var(3,:,:,:));
    end
  end
end
% number of dimensions, neglecting array dimensions of size 1. 
%i.e.  128x128x1 will be 2 dimensions.

dim = ndims(U.bx) - numel(find(size(U.bx)==1)); 
ncell=max(size(U.bx));

idx = ncell/U.L;
bxhat=fftn(U.bx);
byhat=fftn(U.by);
bzhat=fftn(U.bz);

Bspec = 0.5*real( bxhat.*conj(bxhat) + byhat.*conj(byhat) + bzhat.*conj(bzhat))/ncell^(2*dim);

clearvars bxhat byhat bzhat

Bsqr = U.bx.*U.bx + U.by.*U.by + U.bz.*U.bz;
Bsqrhat = fftn(Bsqr);
Bsqrspec = 0.5*real(Bsqrhat.*conj(Bsqrhat))/ncell^(2*dim);
clearvars Bsqr Bsqrhat

%calculate B dot nabla B
gb11 = 0.5*idx*(circshift(U.bx,-1,1) - circshift(U.bx,1,1));
gb21 = 0.5*idx*(circshift(U.bx,-1,2) - circshift(U.bx,1,2));
gb31 = 0.5*idx*(circshift(U.bx,-1,3) - circshift(U.bx,1,3));
gb12 = 0.5*idx*(circshift(U.by,-1,1) - circshift(U.by,1,1));
gb22 = 0.5*idx*(circshift(U.by,-1,2) - circshift(U.by,1,2));
gb32 = 0.5*idx*(circshift(U.by,-1,3) - circshift(U.by,1,3));
gb13 = 0.5*idx*(circshift(U.bz,-1,1) - circshift(U.bz,1,1));
gb23 = 0.5*idx*(circshift(U.bz,-1,2) - circshift(U.bz,1,2));
gb33 = 0.5*idx*(circshift(U.bz,-1,3) - circshift(U.bz,1,3));

% B dot Nabla B
bdgbx = U.bx.*gb11 + U.by.*gb21 + U.bz.*gb31; 
bdgby = U.bx.*gb12 + U.by.*gb22 + U.bz.*gb32;
bdgbz = U.bx.*gb13 + U.by.*gb23 + U.bz.*gb33;

Bdgbxhat=fftn(bdgbx);
Bdgbyhat=fftn(bdgby);
Bdgbzhat=fftn(bdgbz);
clearvars bdgbx bdgby bdgbz

Bdgbspec = 0.5*real( Bdgbxhat.*conj(Bdgbxhat) + Bdgbyhat.*conj(Bdgbyhat) + Bdgbzhat.*conj(Bdgbzhat))/ncell^(2*dim);
clearvars Bdgbxhat Bdgbyhat Bdgbzhat  

% Calculate J
jx = gb23 - gb32; 
jy = gb31 - gb13;
jz = gb12 - gb21;

% B dot J
BdJ = U.bx.*jx + U.by.*jy + U.bz.*jz;
BdJhat=fftn(BdJ);
BdJspec = 0.5*real( BdJhat.*conj(BdJhat))/ncell^(2*dim);

clearvars BdJ BdJhat

% B cross J
[BcJx,BcJy, BcJz]= cross(U.bx,U.by,U.bz,jx,jy,jz);
clearvars jx jy jz 

BcJxhat=fftn(BcJx);
BcJyhat=fftn(BcJy);
BcJzhat=fftn(BcJz);
clearvars BcJx BcJy BcJz  

BcJspec = 0.5*real( BcJxhat.*conj(BcJxhat) + BcJyhat.*conj(BcJyhat) + BcJzhat.*conj(BcJzhat))/ncell^(2*dim);
clearvars BcJxhat BcJyhat BcJzhat  

[~, Bspec_arr]    =calc_spectrum(Bspec,U.L,compensated);
[~, Bsqrspec_arr] =calc_spectrum(Bsqrspec,U.L,compensated);
[~, BcJspec_arr]  =calc_spectrum(BcJspec,U.L,compensated);
[~, BdJspec_arr]  =calc_spectrum(BdJspec,U.L,compensated);
[range, BdgBspec_arr] =calc_spectrum(Bdgbspec,U.L,compensated);

l=length(range);
specout = fopen(sprintf('%s/B_spec_%04.0f.dat',basename,U.time),'w');
for i = 1:l
    fprintf(specout,'%e %e ',range(i),Bspec_arr(i));
    fprintf(specout,'%e %e ',Bsqrspec_arr(i),BcJspec_arr(i));
    fprintf(specout,'%e %e ',BdJspec_arr(i),BdgBspec_arr(i));
    fprintf(specout,'\n');
end
fclose(specout);

function [cx,cy,cz] = cross(ax,ay,az,bx,by,bz)
    cx = ay.*bz - az.*by;
    cy = az.*bx - ax.*bz;
    cz = ax.*by - ay.*bx;
end

end
