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

L=16000;

basename='/tigress/dstonge/pegasus/huge_run2';
if(nargin ==2)
   basename = basename_in;
end    

I=sqrt(-1);


if(nargin ==0)

    % Read in velocity field   VTK
    fname = '/tigress/dstonge/pegasus/huge_run2/combined/combined.0034.fld.vtk'; % filename

    % open file and initialize grid
    mom_ary = [3 3];
    [Grid,status] = init_grid(fname,mom_ary);

    L = Grid.x1max - Grid.x1min;

    for varid=1:length(mom_ary)
        [U.time,var,name,status] = readvtk(Grid,fname,varid);
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


bsqr  = U.bx.^2 + U.by.^2 + U.bz.^2;
bmag = sqrt(bsqr);    
bslices = zeros(3,ncell,ncell);

bslices(1,:,:)=bmag(1,:,:);
bslices(2,:,:)=bmag(:,1,:);
bslices(3,:,:)=bmag(:,:,1);

for i = 1:3
    fslice = fopen(sprintf('%s/B_slice%d_%04.0f.dat',basename,i,U.time),'w');
    for j = 1:ncell
        for k=1:ncell
            fprintf(fslice,'%e %e %e \n',j,k, bslices(i,j,k));
        end
        fprintf(fslice,'\n');
    end
    fclose(fslice);
end

    
end
