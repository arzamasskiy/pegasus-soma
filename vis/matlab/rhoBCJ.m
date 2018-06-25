function rhoBCJ(U,basename_in)
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

ncell=max(size(U.bx));
idx = ncell/U.L;

bmag = sqrt(U.bx.^2 + U.by.^2 + U.bz.^2);
nx = U.bx./bmag;
ny = U.by./bmag;
nz = U.bz./bmag;

%calculate bb dot nable u
gb21 = 0.5*idx*(circshift(U.bx,-1,2) - circshift(U.bx,1,2));
gb31 = 0.5*idx*(circshift(U.bx,-1,3) - circshift(U.bx,1,3));
gb12 = 0.5*idx*(circshift(U.by,-1,1) - circshift(U.by,1,1));
gb32 = 0.5*idx*(circshift(U.by,-1,3) - circshift(U.by,1,3));
gb13 = 0.5*idx*(circshift(U.bz,-1,1) - circshift(U.bz,1,1));
gb23 = 0.5*idx*(circshift(U.bz,-1,2) - circshift(U.bz,1,2));

jx = gb23 - gb32; 
jy = gb31 - gb13;
jz = gb12 - gb21;

clearvars gb12 gb13 gb21 gb23 gb31 gb32 

% B cross J
[BcJx,BcJy, BcJz]= cross(U.bx,U.by,U.bz,jx,jy,jz);
BcJ_mag = sqrt(dot(BcJx,BcJy,BcJz,BcJx,BcJy,BcJz));

clearvars jx jy jz 
clearvars BcJx BcJy BcJz

P = (U.pxx + U.pyy + U.pzz)/3.0;
P_par = nx.*(U.pxx.*nx + U.pxy.*ny + U.pxz.*nz) ...
      + ny.*(U.pxy.*nx + U.pyy.*ny + U.pyz.*nz) ...
      + nz.*(U.pxz.*nx + U.pyz.*ny + U.pzz.*nz);
P_perp = 1.5*P- 0.5*P_par;
rho = sqrt(2.0*P_perp./U.n)./bmag;
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'k_bcj_rho',U.time),sqrt(BcJ_mag(:)./bmag(:)),rho(:),true,true);

clearvars P P_par P_perp bmag rho
clearvars BcJ_mag
    
end
