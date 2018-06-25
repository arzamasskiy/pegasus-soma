function wavenumbers(U,basename_in)
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

fp1 = fopen([basename '/wavenumbers.dat'],write);

if(U.time == 0)
  index=1;
  fprintf(fp1,'[%d]time    ',index) ; index = index+1;
  fprintf(fp1,'[%d]B2   ',index) ; index = index+1;
  fprintf(fp1,'[%d]B4   ',index) ; index = index+1;
  fprintf(fp1,'[%d]U2   ',index) ; index = index+1;
  fprintf(fp1,'[%d]B2U2 ',index) ; index = index+1;
  fprintf(fp1,'[%d]k_par   ',index) ; index = index+1;
  fprintf(fp1,'[%d]k_rms   ',index) ; index = index+1;
  fprintf(fp1,'[%d]k_bcj   ',index) ; index = index+1;
  fprintf(fp1,'[%d]k_bdj   ',index) ; index = index+1;
  fprintf(fp1,'[%d]K_rms   ',index) ; index = index+1;
  fprintf(fp1,'[%d]KB_rms  ',index) ; index = index+1;
  fprintf(fp1,'[%d]Bdnu_sq ',index) ; index = index+1;
  fprintf(fp1,'[%d]bdnu_sq ',index) ; index = index+1;
  fprintf(fp1,'[%d]BBnu    ',index) ; index = index+1;
  fprintf(fp1,'[%d]BBnu_sq ',index) ; index = index+1;
  fprintf(fp1,'[%d]bbnu    ',index) ; index = index+1;
  fprintf(fp1,'[%d]bbnu_sq ',index) ; index = index+1;
  fprintf(fp1,'[%d]nu:nu   ',index) ; index = index+1;
  fprintf(fp1,'[%d]Upar_sq ',index) ; index = index+1;
  fprintf(fp1,'[%d]BdU_sq  ',index) ; index = index+1;
  fprintf(fp1,'[%d]udJ     ',index) ; index = index+1;
  fprintf(fp1,'[%d]udBcJ   ',index) ; index = index+1;
  fprintf(fp1,'[%d]BdJ     ',index) ; index = index+1;
  fprintf(fp1,'[%d]Jrms    ',index) ; %index = index+1;    
  fprintf(fp1,'\n');
end

fprintf(fp1,'%e ',U.time);

%dk = 2 * pi / L;
%kmax = (ncell/2 -1)*dk;
%kmin = - (ncell/2 )*dk;
%[kx,ky,kz] = ndgrid(kmin:dk:kmax);
%kx = ifftshift(kx);
%ky = ifftshift(ky);
%kz = ifftshift(kz);

%ksqr = kx.^2 + ky.^2 + kz.^2;
%k=sqrt(ksqr);

%clearvars kx ky kz ksqr;

% Get rid of the asymmetric Fourier term. Ensures reality.
%chop = ones(ncell,ncell,ncell);
%chop(1,:,:) = zeros(ncell);
%chop(:,1,:) = zeros(ncell);
%chop(:,:,1) = zeros(ncell);
%chop = ifftshift(chop);



usqr = U.vx.^2 + U.vy.^2 + U.vz.^2;
bsqr = U.bx.^2 + U.by.^2 + U.bz.^2;

bmag = sqrt(bsqr);    
bfrth = bsqr.^2;
BU2 = bsqr.*usqr;

B2m = mean(bsqr(:));
B4m = mean(bfrth(:));
U2m = mean(usqr(:));
BU2m = mean(BU2(:));

fprintf(fp1,'%e %e %e %e ',B2m,B4m,U2m,BU2m);

%unit vectors for K_rms (b and n are equivalent)
nx = U.bx./bmag;
ny = U.by./bmag;
nz = U.bz./bmag;

clearvars  bfrth BU2 usqr


%calculate the dyad NABLA B
%bxhat = fftn(U.bx);
%byhat = fftn(U.by);
%bzhat = fftn(U.bz);
%[gb11_h,gb12_h,gb13_h,gb21_h,gb22_h,gb23_h,gb31_h,gb32_h,gb33_h] = dyad(I*chop.*kx,I*chop.*ky,I*chop.*kz,bxhat,byhat,bzhat);

%gb11 = ifftn(gb11_h);
%gb12 = ifftn(gb12_h); 
%gb13 = ifftn(gb13_h);
%gb21 = ifftn(gb21_h);
%gb22 = ifftn(gb22_h);
%gb23 = ifftn(gb23_h);
%gb31 = ifftn(gb31_h);
%gb32 = ifftn(gb32_h);
%gb33 = ifftn(gb33_h);

%clearvars gb11_h gb12_h gb13_h gb21_h gb22_h gb23_h gb31_h gb32_h gb33_h

%calculate bb dot nable u
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
bdgb_sq = abs(dot(bdgbx,bdgby,bdgbz,bdgbx,bdgby,bdgbz));

% tensor contraction of Nabla B with itself
gb_ip = abs( gb11.*gb11 + gb12.*gb12 + gb13.*gb13 ...
           + gb21.*gb21 + gb22.*gb22 + gb23.*gb23 ...
           + gb31.*gb31 + gb32.*gb32 + gb33.*gb33);


k_par = sqrt(mean(bdgb_sq(:))/ B4m);
k_rms = sqrt(mean(gb_ip(:))  / B2m);
       
fprintf(fp1,'%e %e ',k_par,k_rms);

%calculate K
Tdn = nx.*bdgbx + ny.*bdgby + nz.*bdgbz;
Kx = (bdgbx - Tdn.*nx)./bsqr;
Ky = (bdgby - Tdn.*ny)./bsqr;
Kz = (bdgbz - Tdn.*nz)./bsqr;


Ksqr = Kx.*Kx + Ky.*Ky + Kz.*Kz;

K_rms = sqrt(mean(Ksqr(:)));
KB_rms = sqrt(mean(bsqr(:).*Ksqr(:)));

clearvars bdgnx bdgny bdgnz bgdn_sq gb_ip
clearvars Tdn Kx Ky Kz Ksqr bsqr

% Calculate J
%[jx_hat,jy_hat,jz_hat] = cross(I*chop.*kx,I*chop.*ky,I*chop.*kz,bxhat,byhat,bzhat);

%clearvars bxhat byhat bzhat

%jx = ifftn(ifftshift(jx_hat));
%jy = ifftn(ifftshift(jy_hat));
%jz = ifftn(ifftshift(jz_hat));

%clearvars jx_hat jy_hat jz_hat

jx = gb23 - gb32; 
jy = gb31 - gb13;
jz = gb12 - gb21;

clearvars gb11 gb12 gb13 gb21 gb22 gb23 gb31 gb32 gb33

% B dot J
BdJ = dot(U.bx,U.by,U.bz,jx,jy,jz);
BdJ_sq = abs(BdJ.*BdJ);

% B cross J
[BcJx,BcJy, BcJz]= cross(U.bx,U.by,U.bz,jx,jy,jz);
BcJ_sq = abs(dot(BcJx,BcJy,BcJz,BcJx,BcJy,BcJz));

k_bcj = sqrt(mean(BcJ_sq(:)) / B4m);
k_bdj = sqrt(mean(BdJ_sq(:)) / B4m);

ba_BdJ = mean(BdJ(:));

udJ = U.vx.*jx + U.vy.*jy  + U.vz.*jz;
udBcJ = U.vx.*BcJx + U.vy.*BcJy + U.vz.*BcJz;
ba_udJ = mean(udJ(:));
ba_udBcJ = mean(udBcJ(:));

jsqr = jx.*jx + jy.*jy + jz.*jz;
jrms = sqrt(mean(jsqr(:)));

fprintf(fp1,'%e %e ',k_bcj,k_bdj);

clearvars jx jy jz BdJ BcJx BcJy BcJz BdJ_sq BcJ_sq udJ udBcJ jsqr

%calculate the dyad NABLA b
%nxhat = fftn(nx);
%nyhat = fftn(ny);
%nzhat = fftn(nz);

%[gn11_h,gn12_h,gn13_h,gn21_h,gn22_h,gn23_h,gn31_h,gn32_h,gn33_h] = dyad(I*chop.*kx,I*chop.*ky,I*chop.*kz,nxhat,nyhat,nzhat);

%clearvars nxhat nyhat nzhat

%gn11 = ifftn(gn11_h);
%gn12 = ifftn(gn12_h);
%gn13 = ifftn(gn13_h);
%gn21 = ifftn(gn21_h);
%gn22 = ifftn(gn22_h);
%gn23 = ifftn(gn23_h);
%gn31 = ifftn(gn31_h);
%gn32 = ifftn(gn32_h);
%gn33 = ifftn(gn33_h);

%clearvars gn11_h gn12_h gn13_h gn21_h gn22_h gn23_h gn31_h gn32_h gn33_h

%gn11 = 0.5*idx*(circshift(nx,-1,1) - circshift(nx,1,1));
%gn21 = 0.5*idx*(circshift(nx,-1,2) - circshift(nx,1,2));
%gn31 = 0.5*idx*(circshift(nx,-1,3) - circshift(nx,1,3));
%gn12 = 0.5*idx*(circshift(ny,-1,1) - circshift(ny,1,1));
%gn22 = 0.5*idx*(circshift(ny,-1,2) - circshift(ny,1,2));
%gn32 = 0.5*idx*(circshift(ny,-1,3) - circshift(ny,1,3));
%gn13 = 0.5*idx*(circshift(nz,-1,1) - circshift(nz,1,1));
%gn23 = 0.5*idx*(circshift(nz,-1,2) - circshift(nz,1,2));
%gn33 = 0.5*idx*(circshift(nz,-1,3) - circshift(nz,1,3));

%trgb = real(gb11 + gb22 + gb33); %sanity check? Should be zero.


% b dot Nabla b (b or n are the magnetic field unit vector)
%ndgnx = nx.*gn11 + ny.*gn21 + nz.*gn31; 
%ndgny = nx.*gn12 + ny.*gn22 + nz.*gn32;
%ndgnz = nx.*gn13 + ny.*gn23 + nz.*gn33;
%ndgn_sq = abs(dot(ndgnx,ndgny,ndgnz,ndgnx,ndgny,ndgnz));


%clearvars gn11 gn12 gn13 gn21 gn22 gn23 gn31 gn32 gn33
%clearvars ndgnx ndgny ndgnz ngdn_sq


fprintf(fp1,'%e %e ',K_rms, KB_rms);


%calculate the dyad NABLA u
%uxhat = fftn(U.vx);
%uyhat = fftn(U.vy);
%uzhat = fftn(U.vz);

%[gu11_h,gu12_h,gu13_h,gu21_h,gu22_h,gu23_h,gu31_h,gu32_h,gu33_h] = dyad(I*chop.*kx,I*chop.*ky,I*chop.*kz,uxhat,uyhat,uzhat);

%clearvars uxhat uyhat uzhat

%gu11 = ifftn(gu11_h);
%gu12 = ifftn(gu12_h); 
%gu13 = ifftn(gu13_h);
%gu21 = ifftn(gu21_h);
%gu22 = ifftn(gu22_h);
%gu23 = ifftn(gu23_h);
%gu31 = ifftn(gu31_h);
%gu32 = ifftn(gu32_h);
%gu33 = ifftn(gu33_h);

%clearvars gu11_h gu12_h gu13_h gu21_h gu22_h gu23_h gu31_h gu32_h gu33_h

gu11 = 0.5*idx*(circshift(U.vx,-1,1) - circshift(U.vx,1,1));
gu21 = 0.5*idx*(circshift(U.vx,-1,2) - circshift(U.vx,1,2));
gu31 = 0.5*idx*(circshift(U.vx,-1,3) - circshift(U.vx,1,3));
gu12 = 0.5*idx*(circshift(U.vy,-1,1) - circshift(U.vy,1,1));
gu22 = 0.5*idx*(circshift(U.vy,-1,2) - circshift(U.vy,1,2));
gu32 = 0.5*idx*(circshift(U.vy,-1,3) - circshift(U.vy,1,3));
gu13 = 0.5*idx*(circshift(U.vz,-1,1) - circshift(U.vz,1,1));
gu23 = 0.5*idx*(circshift(U.vz,-1,2) - circshift(U.vz,1,2));
gu33 = 0.5*idx*(circshift(U.vz,-1,3) - circshift(U.vz,1,3));


% B dot Nabla U
Bdgux = U.bx.*gu11 + U.by.*gu21 + U.bz.*gu31; 
Bdguy = U.bx.*gu12 + U.by.*gu22 + U.bz.*gu32;
Bdguz = U.bx.*gu13 + U.by.*gu23 + U.bz.*gu33;
Bdgu_sq = abs(dot(Bdgux,Bdguy,Bdguz,Bdgux,Bdguy,Bdguz));

% b dot Nabla U
bdgux = nx.*gu11 + ny.*gu21 + nz.*gu31; 
bdguy = nx.*gu12 + ny.*gu22 + nz.*gu32;
bdguz = nx.*gu13 + ny.*gu23 + nz.*gu33;
bdgu_sq = abs(dot(bdgux,bdguy,bdguz,bdgux,bdguy,bdguz));

fprintf(fp1,'%e %e ',mean(Bdgu_sq(:)),mean(bdgu_sq(:)));

% BB ddot Nabla U
BBdnu  = U.bx.*Bdgux + U.by.*Bdguy + U.bz.*Bdguz; 
BBdnusq = BBdnu.^2; 

% bb ddot Nabla U
bbdnu  = nx.*bdgux + ny.*bdguy + nz.*bdguz; 
bbdnusq = bbdnu.^2; 

fprintf(fp1,'%e %e %e %e ',mean(BBdnu(:)),mean(BBdnusq(:)),mean(bbdnu(:)),mean(bbdnusq(:)));

clearvars Bdgux Bdguy Bdguz bdgux bdguy bdguz BBdnu

% tensor contraction of Nabla U with itself
gu_ip = abs( gu11.*gu11 + gu12.*gu12 + gu13.*gu13 ...
           + gu21.*gu21 + gu22.*gu22 + gu23.*gu23 ...
           + gu31.*gu31 + gu32.*gu32 + gu33.*gu33);

clearvars gu11 gu12 gu13 gu21 gu22 gu23 gu31 gu32 gu33


fprintf(fp1,'%e ',mean(gu_ip(:)));

Upar2  =(nx.*U.vx + ny.*U.vy + nz.*U.vz).^2;
BUpar2 =( U.bx.*U.vx + U.by.*U.vy + U.bz.*U.vz).^2;


fprintf(fp1,'%e %e ',mean(Upar2(:)),mean(BUpar2(:)));


fprintf(fp1,'%e %e %e %e ',ba_udJ, ba_udBcJ,ba_BdJ, jrms);

clearvars gu_ip Bdgu_sq bdgu_sq
clearvars nx ny nz

fprintf(fp1,'\n');
fclose(fp1);



    
end
