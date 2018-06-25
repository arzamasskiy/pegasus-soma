function PDF(U,basename)
%clear all;  

if(nargin ==0)
% Read in velocity moment VTK
basename='/tigress/dstonge/pegasus/sat_b4_proto_3/';
fname  = [basename 'combined/combined.0081.mom.vtk']; % filename

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
fname = [basename 'combined/combined.0100.fld.vtk']; % filename

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

write = 'a';
if(U.time == 0) 
  write = 'w';
end

%fp1 = fopen([basename '~/temp.dat'],write);
fp1 = fopen([basename '/thresholds.dat'],write);


if(U.time == 0)
  index=1;
  fprintf(fp1,'[%d]time                 ',index) ; index = index+1;
  fprintf(fp1,'[%d]beta_perp del t      ',index) ; index = index+1;
  fprintf(fp1,'[%d](1/(beta_perp Delt)) ',index) ; index = index+1;
  fprintf(fp1,'[%d]1/beta_perp          ',index) ; index = index+1;
  fprintf(fp1,'[%d]delt^2 B             ',index) ; index = index+1;
  fprintf(fp1,'[%d]1/ delt^2 B          ',index) ; index = index+1;
  fprintf(fp1,'[%d]firehose per         ',index) ; index = index+1;
  fprintf(fp1,'[%d]mirror per           ',index) ; index = index+1;
  fprintf(fp1,'[%d]Lar median           ',index) ; index = index+1;
  fprintf(fp1,'[%d]Lar mean             ',index) ; index = index+1;
  fprintf(fp1,'[%d]Lar mode             ',index) ; %index = index+1;
  fprintf(fp1,'\n');
end

fprintf(fp1,'%e ', U.time);

%dim = ndims(U.bx) - numel(find(size(U.bx)==1)); 
ncell=max(size(U.bx));
dx = U.L/ncell;
idx = 1.0/dx;

%dk = 2 * pi / U.L;
%kmax = (ncell/2 -1)*dk;
%kmin = - (ncell/2 )*dk;
%[kx,ky,kz] = ndgrid(kmin:dk:kmax);
%kx = ifftshift(kx);
%ky = ifftshift(ky);
%kz = ifftshift(kz);

%I = sqrt(-1);
%ksqr = kx.^2 + ky.^2 + kz.^2;
%k=sqrt(ksqr);

% Get rid of the asymmetric Fourier term. Ensures reality.
%chop = ones(ncell,ncell,ncell);
%chop(1,:,:) = zeros(ncell);
%chop(:,1,:) = zeros(ncell);
%chop(:,:,1) = zeros(ncell);
%chop = ifftshift(chop);


Bsqr = (U.bx).^2 + (U.by).^2 + (U.bz).^2;
Bmag = sqrt(Bsqr);
Brms = sqrt(mean(Bsqr(:)));

%bxhat = fftn(U.bx);
%byhat = fftn(U.by);
%bzhat = fftn(U.bz);
%calculate the dyad NABLA B
%[gb11_h,gb12_h,gb13_h,gb21_h,gb22_h,gb23_h,gb31_h,gb32_h,gb33_h] = dyad(I*chop.*kx,I*chop.*ky,I*chop.*kz,bxhat,byhat,bzhat);

%clearvars bxhat byhat bzhat

%gb11 = ifftn(gb11_h);
%gb12 = ifftn(gb12_h); 
%gb13 = ifftn(gb13_h);
%gb21 = ifftn(gb21_h);
%gb22 = ifftn(gb22_h);
%gb23 = ifftn(gb23_h);
%gb31 = ifftn(gb31_h);
%gb32 = ifftn(gb32_h);
%gb33 = ifftn(gb33_h);

%clearvars gn11_h gb12_h gb13_h gb21_h gb22_h gb23_h gb31_h gb32_h gb33_h

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

%unit vectors for K_rms (b and n are equivalent)
nx = U.bx./Bmag;
ny = U.by./Bmag;
nz = U.bz./Bmag;


%calculate K
Tdn = nx.*bdgbx + ny.*bdgby + nz.*bdgbz;
Kx = (bdgbx - Tdn.*nx)./Bsqr;
Ky = (bdgby - Tdn.*ny)./Bsqr;
Kz = (bdgbz - Tdn.*nz)./Bsqr;


Kmag = sqrt(Kx.*Kx + Ky.*Ky + Kz.*Kz)*U.L/2/pi;

clearvars bdgbx bdgby bdgbz
clearvars Kx Ky Kz


%calculate bb dot nable u
ux_x = 0.5*idx*(circshift(U.vx,-1,1) - circshift(U.vx,1,1));
ux_y = 0.5*idx*(circshift(U.vx,-1,2) - circshift(U.vx,1,2));
ux_z = 0.5*idx*(circshift(U.vx,-1,3) - circshift(U.vx,1,3));

bbux = U.bx.*(U.bx.*ux_x  + U.by.*ux_y + U.bz.*ux_z);
clearvars ux_x ux_y ux_z

uy_x = 0.5*idx*(circshift(U.vy,-1,1) - circshift(U.vy,1,1));
uy_y = 0.5*idx*(circshift(U.vy,-1,2) - circshift(U.vy,1,2));
uy_z = 0.5*idx*(circshift(U.vy,-1,3) - circshift(U.vy,1,3));

bbuy = U.by.*(U.bx.*uy_x  + U.by.*uy_y + U.bz.*uy_z);
clearvars uy_x uy_y uy_z

uz_x = 0.5*idx*(circshift(U.vz,-1,1) - circshift(U.vz,1,1));
uz_y = 0.5*idx*(circshift(U.vz,-1,2) - circshift(U.vz,1,2));
uz_z = 0.5*idx*(circshift(U.vz,-1,3) - circshift(U.vz,1,3));

bbuz = U.bz.*(U.bx.*uz_x  + U.by.*uz_y + U.bz.*uz_z);
clearvars uz_x uz_y uz_z


BBnu = bbux + bbuy + bbuz;

clearvars bbux bbuy bbuz 

bbnu = BBnu./Bsqr;

clearvars BBnu

P = (U.pxx + U.pyy + U.pzz)/3.0;
P_par = nx.*(U.pxx.*nx + U.pxy.*ny + U.pxz.*nz) ...
      + ny.*(U.pxy.*nx + U.pyy.*ny + U.pyz.*nz) ...
      + nz.*(U.pxz.*nx + U.pyz.*ny + U.pzz.*nz);
P_perp = 1.5*P- 0.5*P_par;

vperp = sqrt(2.0*P_perp./U.n);

beta_perp = 2 * P_perp ./ Bsqr;
beta_par = 2 * P_par ./ Bsqr;
beta = 2 * P./Bsqr;


Delt = P_perp./P_par - 1.0;

lar = sqrt(2*P_perp./(U.n.*Bsqr));


ba_thres1 = mean(Delt(:));
Delt_bp = Delt.*beta_perp;
ba_beta_perp = 1.0./(mean(1.0./beta_perp(:)));
ba_thres2 = mean(Delt_bp(:));

mg= Delt.^2.*Bmag;
ba_m_growth1 = mean(mg(:));
ba_m_growth2 = 1.0./mean(1.0./mg(:));

fprintf(fp1,'%e %e %e ',ba_thres1, ba_thres2, ba_beta_perp);
fprintf(fp1,'%e %e ',ba_m_growth1, ba_m_growth2);

%clearvars Bsqr Delt_bp mg

clearvars P P_par P_perp
clearvars nx ny nz

% Create our histograms
[~, ~, ~] = mkdir(sprintf('%s/PDFs',basename));


ind=Delt > (1.0./beta_perp);
nm = length(find(ind(:)));
ind=Delt < -(2./beta_par);
nf = length(find(ind(:)));
tot = length(beta_perp(:));

fprintf(fp1,'%e %e ',nf/tot, nm/tot);

mean_beta=mean(log10(beta(:)));
std_beta=std(log10(beta(:)));
ind= abs(log10(beta)-mean_beta) < std_beta;

print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'DeltxBeta',U.time),Delt(ind).*beta(ind),false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'bbnu_deltxbeta',U.time),bbnu(ind),Delt(ind).*beta(ind),false,false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'bbnu_rho',U.time),bbnu(:),vperp(:)./Bmag(:),false,true);


print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'B',U.time),Bmag(:),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'beta',U.time),beta(:),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'B_Brms',U.time),Bmag(:)/Brms,true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag',U.time),Kmag(:),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'delt',U.time),Delt(:),false);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'ROS',U.time),bbnu(:),false);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'lar',U.time),lar(:),true);

[hist,edges]=histcounts(lar(:));

l=length(hist(:));
pdfbin=hist.*(edges(2:(l+1)) - edges(1:l));
range = 0.5*(edges(2:(l+1)) + edges(1:l));
ave = sum(pdfbin(:).*range(:))/sum(pdfbin(:));
[~,ind] = max(pdfbin(:));
mode = range(ind);
median=0;
s=0;
pdfbin  = pdfbin/sum(pdfbin(:));
for i=1:l
    s = s + pdfbin(i);
    if(s>0.5)
       median = range(i);
       break;
    end
end


fprintf(fp1,'%e %e %e ',median, ave, mode);

print_2D_PDF(sprintf('%s/mirror_fire_%04.0f.dat',basename,U.time),beta_par(:),Delt(:),true,false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Delt_ROS',U.time),Delt(:),bbnu(:),false,false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_B_Brms',U.time),Kmag(:),Bmag(:)/Brms,true,true);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_ROS',U.time),Kmag(:),bbnu(:),true,false);



ind=(Delt > 0);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_P',U.time),Kmag(ind),true);

ind=(Delt < 0);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_N',U.time),Kmag(ind),true);

ind=Delt > (1.0./beta_perp);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_M',U.time),Kmag(ind),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'ROS_M',U.time),bbnu(ind),false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_B_Brms_M',U.time),Kmag(ind),Bmag(ind)/Brms,true,true);

ind=Delt < -(2./beta_par);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_F',U.time),Kmag(ind),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'ROS_F',U.time),bbnu(ind),false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_B_Brms_F',U.time),Kmag(ind),Bmag(ind)/Brms,true,true);

ind=(Delt >= -(2./beta_par)) & (Delt <= (1.0./beta_perp));
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_S',U.time),Kmag(ind),true);
print_1D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'ROS_S',U.time),bbnu(ind),false);
print_2D_PDF(sprintf('%s/PDFs/PDF_%s_%04.0f.dat',basename,'Kmag_B_Brms_S',U.time),Kmag(ind),Bmag(ind)/Brms,true,true);

fprintf(fp1,'\n');

fclose(fp1);

return;

end

