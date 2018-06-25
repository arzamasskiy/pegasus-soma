function V_spectrum(U,basename_in)
%
% Plots magnetic field spectrum. 
% Tested for 3 dimensions. 
% Should work for any dimension
%
% Also will plot relevant wave-numbers 
% given in 'Simulations of the Small-Scale Turbulent Dynamo',
% Schekochihin et. al. 2004.
%

compensated=true; % try to compensate for grid cut-offs
basename='/tigress/dstonge/pegasus/prod_run3_FP5/';


if(nargin ==2)
   basename = basename_in;
end    

if(nargin ==0)
  % Read in velocity field   VTK
  fname = '/tigress/dstonge/pegasus/prod_run3_FP5/combined/combined.0000.mom.vtk'; % filename
  %fname = '/tigress/dstonge/pegasus/run7_alt2/combined/dynamo.0050.fld.vtk'; % filename
  % open file and initialize grid
  mom_ary = [1 3 1 1 1 1 1 1];
  [Grid,~] = init_grid(fname,mom_ary);

  U.L = Grid.x1max - Grid.x1min;

  for varid=1:length(mom_ary)
    [time,var,~,~] = readvtk(Grid,fname,varid);
    if (varid==1)         % read density
      U.n = squeeze(var);
    elseif (varid==2)    % read momentum density
      U.vx = squeeze(var(1,:,:,:))./U.n;
      U.vy = squeeze(var(2,:,:,:))./U.n;
      U.vz = squeeze(var(3,:,:,:))./U.n;
    elseif (varid==3)    % read pressure tensor
      U.pxx= squeeze(var) - U.n.*U.vx.*U.vx;
    elseif (varid==4)    % read pressure tensor
      U.pxy= squeeze(var) - U.n.*U.vx.*U.vy;
    elseif (varid==5)    % read pressure tensor
      U.pxz= squeeze(var) - U.n.*U.vx.*U.vz;
    elseif (varid==6)    % read pressure tensor
      U.pyy= squeeze(var) - U.n.*U.vy.*U.vy;
    elseif (varid==7)    % read pressure tensor
      U.pyz= squeeze(var) - U.n.*U.vy.*U.vz;
    elseif (varid==8)    % read pressure tensor
      U.pzz= squeeze(var) - U.n.*U.vz.*U.vz;
    end
    U.time = time;
  end
end
% number of dimensions, neglecting array dimensions of size 1. 
%i.e.  128x128x1 will be 2 dimensions.


dim = ndims(U.vx) - numel(find(size(U.vx)==1));
ncell=max(size(U.vx));

vxhat = fftn(U.vx);
vyhat = fftn(U.vy);
vzhat = fftn(U.vz);

Vspec = 0.5*real( vxhat.*conj(vxhat) + vyhat.*conj(vyhat) + vzhat.*conj(vzhat))/ncell^(2*dim);

clearvars vxhat vyhat vzhat;

[range, Vspec_arr] = calc_spectrum(Vspec, U.L,compensated);


l=length(range);

loglog(range,Vspec_arr(1:l));
drawnow

specout = fopen(sprintf('%s/V_spec_%04.0f.dat',basename,U.time),'w');
for i = 1:l
    fprintf(specout,'%e %e\n',range(i),Vspec_arr(i));
end
fclose(specout);
    
end
