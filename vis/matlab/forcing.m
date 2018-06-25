function forcing(U,basename_in)
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

write='a';
if(U.time == 0) 
  write='W';
end

fp1 = fopen([basename '/forcing.dat'],write);

if(U.time == 0)
  index=1;
  fprintf(fp1,'[%d]time    ',index) ; index = index+1;
  fprintf(fp1,'[%d]divF   ',index) ; index = index+1;
  fprintf(fp1,'[%d]FdcF    ',index) ; index = index+1;
  fprintf(fp1,'[%d]Frms    ',index) ; index = index+1;
  fprintf(fp1,'[%d]Hrms    ',index) ; %index = index+1;
  fprintf(fp1,'\n');
end

fprintf(fp1,'%e ',U.time);


fsqr = U.fx.^2 + U.fy.^2 + U.fz.^2;

frms = sqrt(mean(fsqr(:)));

%calculate bb dot nable u
gf11 = 0.5*idx*(circshift(U.fx,-1,1) - circshift(U.fx,1,1));
gf21 = 0.5*idx*(circshift(U.fx,-1,2) - circshift(U.fx,1,2));
gf31 = 0.5*idx*(circshift(U.fx,-1,3) - circshift(U.fx,1,3));
gf12 = 0.5*idx*(circshift(U.fy,-1,1) - circshift(U.fy,1,1));
gf22 = 0.5*idx*(circshift(U.fy,-1,2) - circshift(U.fy,1,2));
gf32 = 0.5*idx*(circshift(U.fy,-1,3) - circshift(U.fy,1,3));
gf13 = 0.5*idx*(circshift(U.fz,-1,1) - circshift(U.fz,1,1));
gf23 = 0.5*idx*(circshift(U.fz,-1,2) - circshift(U.fz,1,2));
gf33 = 0.5*idx*(circshift(U.fz,-1,3) - circshift(U.fz,1,3));

hx = gf23 - gf32; 
hy = gf31 - gf13;
hz = gf12 - gf21;

hsqr = hx.*hx + hy.*hy + hz.*hz;
hrms = sqrt(mean(hsqr(:)));

divF = gf11 + gf22 + gf33;
FdcF = U.fx.*hx +  U.fy.*hy +   U.fz.*hz;

ba_divF = mean(divF(:));
ba_FdcF = mean(FdcF(:));

clearvars gf11 gf12 gf13 gf21 gf22 gf23 gf31 gf32 gf33

fprintf(fp1,'%e %e %e %e ',ba_divF, ba_FdcF, frms, hrms);

fprintf(fp1,'\n');
fclose(fp1);

end
