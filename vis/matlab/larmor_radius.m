function larmor_radius(U,basename)
%
% Finds teh cell-averaged Larmor radius.
% Should be the most useful measure of larmor radius vs. grid size.
%

% Read in velocity moment and field VTK

if(nargin ==0)
  fname_mom = '/tigress/dstonge/pegasus/saturation5/combined/combined.0000.mom.vtk'; % filename
  fname_fld = '/tigress/dstonge/pegasus/saturation5/combined/combined.0000.fld.vtk'; % filename
  %fname_mom = '/tigress/dstonge/pegasus/prod_run3/combined/combined.0000.mom.vtk'; % filename
  %fname_fld = '/tigress/dstonge/pegasus/prod_run3/combined/combined.0000.fld.vtk'; % filename

  % open file and initialize grid
  mom_ary = [3 1 3 1 1 1 1 1 1];
  n=1;
  [Grid,~] = init_grid(fname_mom,mom_ary);

  U.L = Grid.x1max - Grid.x1min;

  for varid=1:length(mom_ary)
    [U.time,var,~,~] = readvtk(Grid,fname_mom,varid);
    if (varid==(1+n))         % read density
      U.n = squeeze(var);
    elseif (varid==(2+n))    % read momentum density
      U.vx = squeeze(var(1,:,:,:))./U.n;
      U.vy = squeeze(var(2,:,:,:))./U.n;
      U.vz = squeeze(var(3,:,:,:))./U.n;
    elseif (varid==3+n)    % read pressure tensor
      U.pxx= squeeze(var) - U.n.*U.vx.*U.vx;
    elseif (varid==4+n)    % read pressure tensor
      U.pxy= squeeze(var) - U.n.*U.vx.*U.vy;
    elseif (varid==5+n)    % read pressure tensor
      U.pxz= squeeze(var) - U.n.*U.vx.*U.vz;
    elseif (varid==6+n)    % read pressure tensor
      U.pyy= squeeze(var) - U.n.*U.vy.*U.vy;
    elseif (varid==7+n)    % read pressure tensor
      U.pyz= squeeze(var) - U.n.*U.vy.*U.vz;
    elseif (varid==8+n)    % read pressure tensor
      U.pzz= squeeze(var) - U.n.*U.vz.*U.vz;
    end
  end

  % open file and initialize grid
  mom_ary = [3 3];
  [Grid,~] = init_grid(fname_fld,mom_ary);

  for varid=1:length(mom_ary)
    [U.time,var,~,~] = readvtk(Grid,fname_fld,varid);
    if (varid==1)        % read cell-centered magnetic field
      U.bx = squeeze(var(1,:,:,:));
      U.by = squeeze(var(2,:,:,:));
      U.bz = squeeze(var(3,:,:,:));
    elseif (varid==2)    % read cell-centered electric field
      %  U.ex = squeeze(var(1,:,:,:));
      %  U.ey = squeeze(var(2,:,:,:));
      %  U.ez = squeeze(var(3,:,:,:));
    end
  end
  clear var;
end

write = 'a';
if(U.time == 0) 
  write = 'w';
end

fp1 = fopen([basename '/larmor.dat'],write);

if(U.time == 0)
  index=1;
  fprintf(fp1,'[%d]time         ',index) ; index = index+1;
  fprintf(fp1,'[%d]rho_box      ',index) ; index = index+1;
  fprintf(fp1,'[%d]1/(1/rho_box)',index) ; index = index+1;
  fprintf(fp1,'[%d]T            ',index) ; index = index+1;
  fprintf(fp1,'[%d]v_par        ',index) ; index = index+1;
  fprintf(fp1,'[%d]v_perp       ',index) ; index = index+1;
  fprintf(fp1,'[%d]delt         ',index) ; index = index+1;
  fprintf(fp1,'[%d]deltrms      ',index) ; index = index+1;
  fprintf(fp1,'[%d]B^2 delt     ',index) ; index = index+1;
  fprintf(fp1,'[%d]Bmag         ',index) ; index = index+1;
  fprintf(fp1,'[%d]Bmag_p       ',index) ; index = index+1;
  fprintf(fp1,'[%d]Brms         ',index) ; index = index+1;
  fprintf(fp1,'[%d]mean_lars    ',index) ; index = index+1;
  fprintf(fp1,'[%d]in_out       ',index) ; %index = index+1;
  fprintf(fp1,'\n');
end

fprintf(fp1,'%e ', U.time);

L = U.L;
Bsqr = (U.bx).^2 + (U.by).^2 + (U.bz).^2;
Bmag = sqrt(Bsqr);
bxhat = (U.bx)./Bmag;
byhat = (U.by)./Bmag;
bzhat = (U.bz)./Bmag;

mean_bmag = mean(Bmag(:));
mean_bmagn = mean(Bmag(:).*U.n(:));
mean_bsqr = sqrt(mean(Bsqr(:)));
  

P = (U.pxx + U.pyy + U.pzz)/3.0;
P_par = bxhat.*(U.pxx.*bxhat + U.pxy.*byhat + U.pxz.*bzhat) ...
      + byhat.*(U.pxy.*bxhat + U.pyy.*byhat + U.pyz.*bzhat) ...
      + bzhat.*(U.pxz.*bxhat + U.pyz.*byhat + U.pzz.*bzhat);
P_perp = 1.5*P - 0.5*P_par;

Delt = P_perp./P_par - 1.0;
ba_delt = mean(Delt(:));
ba_deltrms = sqrt(mean(Delt(:).^2));
ba_bsqr_delt = mean(Bsqr(:).*Delt(:));

mean_T  =  mean(P(:)./U.n(:));

T_par = P_par./U.n;
T_perp = P_perp./U.n;
clearvars P P_par bxhat byhat bzhat Delt;
clearvars Bsqr;


v_par = sqrt(2.0*T_par);
v_perp = sqrt(2.0*T_perp);
mean_vpar  = mean(v_par(:));
mean_vperp = mean(v_perp(:));
lar = v_perp./Bmag;
ilar = 1./lar;



clearvars P_perp T_perp v_perp Bmag;

mean_lar=mean(lar(:));
mean_ilar=mean(ilar(:));

fprintf(fp1,'%e %e ', mean_lar, mean_ilar);
fprintf(fp1,'%e %e %e ', mean_T,mean_vpar,mean_vperp);
fprintf(fp1,'%e %e %e ', ba_delt, ba_deltrms,ba_bsqr_delt);
fprintf(fp1,'%e %e ', mean_bmag,mean_bmagn,mean_bsqr);


lar_array=lar(:);
ind=find(lar_array < L);
count = length(ind(:));
mean_lars = mean(lar_array(ind));
in_out=count/length(lar(:));

fprintf(fp1,'%e %e ', mean_lars, in_out);
fprintf(fp1,'\n');
fclose(fp1);

end
