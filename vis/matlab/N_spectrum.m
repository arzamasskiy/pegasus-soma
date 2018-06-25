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
basename='/tigress/dstonge/pegasus/prod_run3';


I=sqrt(-1);

if(nargin ==2)
   basename = basename_in;
end    

if(nargin ==0)
  % Read in velocity field   VTK
  fname = '/tigress/dstonge/pegasus/prod_run3/combined/combined.0050.mom.vtk'; % filename
  %fname = '/tigress/dstonge/pegasus/run7_alt2/combined/dynamo.0050.fld.vtk'; % filename

  % open file and initialize grid
  mom_ary = [1 3 1 1 1 1 1 1];
  [Grid,status] = init_grid(fname,mom_ary);
  U.L = Grid.x1max - Grid.x1min;
  [U.time,var,name,status] = readvtk(Grid,fname,1); %read density
  U.n = squeeze(var);
end
% number of dimensions, neglecting array dimensions of size 1. 
%i.e.  128x128x1 will be 2 dimensions.




dim = ndims(U.n) - numel(find(size(U.n)==1))
ncell=max(size(U.n));

nhat = fftn(U.n);


Nspec = 0.5*real( nhat.*conj(nhat))/ncell^(2*dim);

clearvars nhat;

[range Nspec_arr] = calc_spectrum(nhat,U.L,compensated);
l=length(range);

specout = fopen(sprintf('%s/N_spec_%04.0f.dat',basename,U.time),'w');
for i = 1:l
    fprintf(specout,'%e %e\n',range(i),Nspec_arr(i));
end
fclose(specout);
    
end
