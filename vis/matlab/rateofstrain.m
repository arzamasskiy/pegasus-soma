function rateofstrain(U,basename)
%clear all;  

compensated = false; 

if(nargin ==0)
  % Read in velocity moment VTK
  basename='/tigress/dstonge/pegasus/prod_run3/combined';
  fname = '/tigress/dstonge/pegasus/prod_run3/combined/combined.0014.mom.vtk'; % filename

  % open file and initialize grid
  mom_ary = [1 3 1 1 1 1 1 1];
  [Grid,~] = init_grid(fname,mom_ary);

  U.L = Grid.x1max - Grid.x1min;
  for varid=1:length(mom_ary)
       [U.time,var,~,~] = readvtk(Grid,fname,varid);
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
  end

  % Read in velocity field   VTK
  fname = '/tigress/dstonge/pegasus/prod_run3/combined/combined.0028.fld.vtk'; % filename

  % open file and initialize grid
  mom_ary = [3 3];
  [Grid,~] = init_grid(fname,mom_ary);

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


write='a';
if(U.time == 0) 
    write='W';
end

fp1 = fopen([basename '/ROS.dat'],write);

if(U.time == 0)
  index=1;
  fprintf(fp1,'[%d]time      ',index) ; index = index+1;
  fprintf(fp1,'[%d]ROS       ',index) ; index = index+1;
  fprintf(fp1,'[%d]ROS^2     ',index) ; index = index+1;
  fprintf(fp1,'[%d]divU      ',index) ; index = index+1;
  fprintf(fp1,'[%d]divUBB    ',index) ; index = index+1;
  fprintf(fp1,'[%d]helicity  ',index) ; index = index+1;
  fprintf(fp1,'[%d]stdheli   ',index) ; index = index+1;
  fprintf(fp1,'[%d]1/beta    ',index) ; index = index+1;
  fprintf(fp1,'[%d]Bsqr/beta ',index) ; index = index+1;
  fprintf(fp1,'[%d]ROS beta  ',index) ; index = index+1;
  fprintf(fp1,'[%d]1/(ROSbeta)',index); index = index+1;
  fprintf(fp1,'[%d]Diss      ',index) ; index = index+1;
  fprintf(fp1,'[%d]hDiss     ',index) ; index = index+1;
  fprintf(fp1,'[%d]dBdt/Delt ',index) ; index = index+1;
  fprintf(fp1,'[%d]dBdt beta ',index) ; index = index+1;
  fprintf(fp1,'[%d]<bbnu>    ',index) ; index = index+1;
  fprintf(fp1,'[%d]<bbnu>^2  ',index) ; index = index+1;
  fprintf(fp1,'\n');
end
rest = U.rest;
hrest= U.hrest;

fprintf(fp1,'%e ',U.time);

dim = ndims(U.bx) - numel(find(size(U.bx)==1)); 
ncell=max(size(U.bx));
dx = U.L/ncell;
idx = 1.0/dx;

Bsqr = (U.bx).^2 + (U.by).^2 + (U.bz).^2;
Brms = mean(Bsqr(:));
iBrms=1.0/Brms;
Bmag = sqrt(Bsqr);

ux_x = 0.5*idx*(circshift(U.vx,-1,1) - circshift(U.vx,1,1));
ux_y = 0.5*idx*(circshift(U.vx,-1,2) - circshift(U.vx,1,2));
ux_z = 0.5*idx*(circshift(U.vx,-1,3) - circshift(U.vx,1,3));

BBux = U.bx.*(U.bx.*ux_x  + U.by.*ux_y + U.bz.*ux_z);
divu = ux_x;
hel = U.vy.*ux_z - U.vz.*ux_y;
clearvars ux_x ux_y ux_z

uy_x = 0.5*idx*(circshift(U.vy,-1,1) - circshift(U.vy,1,1));
uy_y = 0.5*idx*(circshift(U.vy,-1,2) - circshift(U.vy,1,2));
uy_z = 0.5*idx*(circshift(U.vy,-1,3) - circshift(U.vy,1,3));

BBuy = U.by.*(U.bx.*uy_x  + U.by.*uy_y + U.bz.*uy_z);
divu = divu + uy_y;
hel = hel - U.vx.*uy_z + U.vz.*uy_x;
clearvars uy_x uy_y uy_z

uz_x = 0.5*idx*(circshift(U.vz,-1,1) - circshift(U.vz,1,1));
uz_y = 0.5*idx*(circshift(U.vz,-1,2) - circshift(U.vz,1,2));
uz_z = 0.5*idx*(circshift(U.vz,-1,3) - circshift(U.vz,1,3));

BBuz = U.bz.*(U.bx.*uz_x  + U.by.*uz_y + U.bz.*uz_z);
divu = divu + uz_z;
hel = hel + U.vx.*uz_y - U.vy.*uz_x;
clearvars uz_x uz_y uz_z


ba_divu = mean(divu(:));
ba_hel = mean(hel(:));
sig_hel = std(hel(:));
divu = 0.5*(divu.*Bsqr);

ba_divuBB = mean(divu(:))*iBrms;
clearvars divu hel

BBnu = BBux + BBuy + BBuz;
ba_bbnu = mean(BBnu(:))*iBrms;
ba_bbnusqr = mean(BBnu(:).^2)*iBrms^2;

bbnu = BBnu./Bsqr;
babbnu_57 = mean(bbnu(:));
babbnu_58 = mean(bbnu(:).^2);

fprintf(fp1,'%e %e %e %e %e %e ',ba_bbnu,ba_bbnusqr,ba_divu,ba_divuBB,ba_hel,sig_hel);

clearvars BBux BBuy BBuz bbnu

nx= U.bx./Bmag;
ny= U.by./Bmag;
nz= U.bz./Bmag;

P = (U.pxx + U.pyy + U.pzz)/3.0;
P_par = nx.*(U.pxx.*nx + U.pxy.*ny + U.pxz.*nz) ...
      + ny.*(U.pxy.*nx + U.pyy.*ny + U.pyz.*nz) ...
      + nz.*(U.pxz.*nx + U.pyz.*ny + U.pzz.*nz);
P_perp = 1.5*P - 0.5*P_par;

Delt = P_perp./P_par - 1.0;
beta = 2*P ./ Bmag.^2;
beta_perp = 2 * P_perp ./ Bsqr;
beta_par = 2 * P_par ./ Bsqr;

Pperphat = fftn(P_perp);
Pparhat = fftn(P_par);

Pperpspec = real( Pperphat.*conj(Pperphat))/ncell^(2*dim);
Pparspec = real( Pparhat.*conj(Pparhat))/ncell^(2*dim);
clearvars Pperphat Pparhat

[~, Pperpspec_arr] = calc_spectrum(Pperpspec, U.L,compensated);
[~, Pperpspec_arr] = calc_spectrum(Pperpspec, U.L,compensated);
[range, Pparspec_arr] = calc_spectrum(Pparspec, U.L,compensated);

l=length(range);

specout = fopen(sprintf('%s/P_spec_%04.0f.dat',basename,U.time),'w');
for i = 1:l
    fprintf(specout,'%e %e %e\n',range(i),Pperpspec_arr(i),Pparspec_arr(i));
end
fclose(specout);

clearvars Pperpspec Pparspec Pperpspec_arr Pparspec_arr range

ba_beta = 1.0./mean(1.0./beta(:));
ba_bsqr_beta = mean(Bsqr(:)./beta(:));


beta_ros = 2*P.*BBnu./Bsqr;
%ba_beta_ros2 = 1.0./(sum(1.0./beta_ros(:))./ncell^dim)
ba_beta_ros1 = mean(beta_ros(:))*iBrms;
ba_beta_ros2 = mean(beta_ros(:)./Bsqr(:));

fprintf(fp1,'%e %e ',ba_beta,ba_bsqr_beta);
fprintf(fp1,'%e %e ',ba_beta_ros1,ba_beta_ros2);


dk = 2 * pi / U.L;
kmax = (ncell/2 -1)*dk;
kmin = - (ncell/2 )*dk;
[kx,ky,kz] = ndgrid(kmin:dk:kmax);
kx = ifftshift(kx);
ky = ifftshift(ky);
kz = ifftshift(kz);

%ksqr = kx.^2 + ky.^2 + kz.^2;
% grid k^2 from Von Neumann stability analysis (what the code actually sees)
ksqr = 2*((1-cos(kx*dx)) + (1-cos(ky*dx))+(1-cos(kz*dx)))/dx^2;

bxhat = fftn(U.bx);
byhat = fftn(U.by);
bzhat = fftn(U.bz);

bsqrh = bxhat.*conj(bxhat) + byhat.*conj(byhat) + bzhat.*conj(bzhat);
bsqr_tot = mean(bsqrh(:));

diss  =     rest*ksqr.*(bxhat.*conj(bxhat) + byhat.*conj(byhat) + bzhat.*conj(bzhat));
dissh = hrest*ksqr.^2.*(bxhat.*conj(bxhat) + byhat.*conj(byhat) + bzhat.*conj(bzhat));


ba_rest = mean(diss(:))/(bsqr_tot);
ba_hrest = mean(dissh(:))/(bsqr_tot);
fprintf(fp1,'%e %e ',ba_rest,ba_hrest);

ind=(Delt < -(2./beta_par)) | (Delt > (1.0./beta_perp));

d2Bx = real(ifftn(ksqr.*bxhat));
d2By = real(ifftn(ksqr.*byhat));
d2Bz = real(ifftn(ksqr.*bzhat));
d4Bx = real(ifftn(ksqr.^2.*bxhat));
d4By = real(ifftn(ksqr.^2.*byhat));
d4Bz = real(ifftn(ksqr.^2.*bzhat));

dBdt = (BBnu ...
       - (rest* (U.bx.*d2Bx+ U.by.*d2By + U.bz.*d2Bz)  ...
       + hrest*(U.bx.*d4Bx + U.by.*d4By + U.bz.*d4Bz)))./Bsqr;
dBdt_Delt = dBdt./Delt;

clearvars d2Bx d2By d2Bz d4Bx d4By d4Bz

dBdtbeta = 2*P.*dBdt./Bsqr;

ba_dBdt_Delt = mean(dBdt_Delt(ind));
ba_dBdtbeta = mean(dBdtbeta(ind));        

fprintf(fp1,'%e %e ',ba_dBdt_Delt,ba_dBdtbeta);

%mean_beta=mean(log10(beta(:)));
%std_beta=std(log10(beta(:)));
%ind= abs(log10(beta)-mean_beta) < std_beta;
%print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'dBdt_deltxbeta',U.time),3*dBdt(ind),Delt(ind).*beta(ind),false,false);
%print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'dBdt_beta',U.time),3*dBdt(:),beta(:),false,true);
%print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'dBdt_delt',U.time),3*dBdt(:),Delt(:),false,false);
%xm = mean(3*dBdt(:));
%ym = mean(Delt(:).*beta(:));
%b = sum((3*dBdt(:) - xm).*(Delt(:).*beta(:)-ym))/sum((3*dBdt(:) - xm).^2);
%b = cov(3*dBdt(:)',Delt(:)'.*beta(:)')/cov(3*dBdt(:)')
%a = ym - b*xm;

clearvars dBdt_Delt dBdtbeta
clearvars bxhat byhat bzhat bsqr_tot;
clearvars kx ky kz ksqr;

BBnu_h = fftn(BBnu/Brms);
bbnu_h = fftn(BBnu./Bsqr);
Delt_h = fftn(Delt);

Vspec = real(BBnu_h.*conj(BBnu_h))/ncell^(dim*2);
Vspec2 = real(bbnu_h.*conj(bbnu_h))/ncell^(dim*2);
Vspec3 = real(Delt_h.*conj(Delt_h))/ncell^(dim*2);

clearvars BBnu_h BBnu_h Delt_h;

fprintf(fp1,'%e %e ',babbnu_57,babbnu_58);

fprintf(fp1,'\n');
fclose(fp1);

[~, Vspec_arr] = calc_spectrum(Vspec,U.L,compensated);
[range, Vspec_arr2] = calc_spectrum(Vspec2,U.L,compensated);
[range, Vspec_arr3] = calc_spectrum(Vspec3,U.L,compensated);

clearvars Vspec Vspec2 Vspec3


specout = fopen(sprintf('%s/ROS_spec_%04.0f.dat',basename,U.time),'w');
for i = 1:length(range(:))
    fprintf(specout,'%e %e %e %e\n',range(i),Vspec_arr(i), Vspec_arr2(i),Vspec_arr3(i));
end
fclose(specout);
    

return;
