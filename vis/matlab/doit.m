clear all; 
t = 0:1000;
cadence_fld=1;
cadence_mom=1;
basename='/tigress/dstonge/pegasus/sat_b4_proto_3/';
diary([basename 'log.txt']);
diary on;
fp1=fopen([basename 'diss.dat'],'r');
diss=fscanf(fp1,'%e %e');
fclose(fp1);
for i = t  
  diary off; diary on;

  clearvars U var time status;
  if(mod(i,cadence_fld)==0 && mod(i,cadence_mom)==0 )
  %if(mod(i,cadence_fld)==0 )
    disp('read fld');
    fname_fld = sprintf([basename '/combined/combined.%04d.fld.vtk'],i/cadence_fld); % filename
    if(exist(fname_fld,'file') ~= 2)
      continue;
    end
    % open file and initialize grid
    mom_ary = [3  ];
    [Grid,status] = init_grid(fname_fld,mom_ary);
    U.L = Grid.x1max - Grid.x1min;
    for varid=1:length(mom_ary)
      [U.time,var,name,status] = readvtk(Grid,fname_fld,varid);
      if (varid==1)    % read cell-centered magnetic field
        U.bx = squeeze(var(1,:,:,:));
        U.by = squeeze(var(2,:,:,:));
        U.bz = squeeze(var(3,:,:,:));
      elseif (varid==2)  % read cell-centered electric field
        U.ex = squeeze(var(1,:,:,:));
        U.ey = squeeze(var(2,:,:,:));
        U.ez = squeeze(var(3,:,:,:));
      elseif (varid==3)  % read cell-centered electric field
        U.fx = squeeze(var(1,:,:,:));
        U.fy = squeeze(var(2,:,:,:));
        U.fz = squeeze(var(3,:,:,:));
      end
    end
    clear var;
    U.rest=diss(1); U.hrest=diss(2);
    disp('Field read')
    diary off; diary on;

     

    if(mod(i,cadence_mom)==0)
      fname_mom = sprintf([basename '/combined/combined.%04d.mom.vtk'],i/cadence_mom); % filename
     
      if(exist(fname_mom,'file') ~= 2)
        continue;
      end
      % open file and initialize grid
     
      disp('read mom');
    
      if(i < 10000)
        mom_ary = [1 3];
        mom_ary = [1 3 1 1 1 1 1 1];
        [Grid,status] = init_grid(fname_mom,mom_ary);
        for varid=1:length(mom_ary)
          [time,var,name,status] = readvtk(Grid,fname_mom,varid);  
          if(abs(time-U.time) > 0.01)
            disp('Time error!');
          end
          if (varid==1)     % read density
            U.n = squeeze(var);
          elseif (varid==2)  % read momentum density
            U.vx = squeeze(var(1,:,:,:))./U.n;
            U.vy = squeeze(var(2,:,:,:))./U.n;
            U.vz = squeeze(var(3,:,:,:))./U.n;
          elseif (varid==3)  % read pressure tensor
            U.pxx= squeeze(var) - U.n.*U.vx.*U.vx;
          elseif (varid==4)  % read pressure tensor
            U.pxy= squeeze(var) - U.n.*U.vx.*U.vy;
          elseif (varid==5)  % read pressure tensor
            U.pxz= squeeze(var) - U.n.*U.vx.*U.vz;
          elseif (varid==6)  % read pressure tensor
            U.pyy= squeeze(var) - U.n.*U.vy.*U.vy;
          elseif (varid==7)  % read pressure tensor
            U.pyz= squeeze(var) - U.n.*U.vy.*U.vz;
          elseif (varid==8)  % read pressure tensor
            U.pzz= squeeze(var) - U.n.*U.vz.*U.vz;
          end
        end
      else
        mom_ary = [3 1 3 1 1 1 1 1 1];
        [Grid,status] = init_grid(fname_mom,mom_ary);
        for varid=1:length(mom_ary)
          [time,var,name,status] = readvtk(Grid,fname_mom,varid);  
          if(abs(time-U.time) > 0.01)
            disp('Time error!');
          end
          if (varid==2)     % read density
            U.n = squeeze(var);
          elseif (varid==3)  % read momentum density
            U.vx = squeeze(var(1,:,:,:))./U.n;
            U.vy = squeeze(var(2,:,:,:))./U.n;
            U.vz = squeeze(var(3,:,:,:))./U.n;
          elseif (varid==4)  % read pressure tensor
            U.pxx= squeeze(var) - U.n.*U.vx.*U.vx;
          elseif (varid==5)  % read pressure tensor
            U.pxy= squeeze(var) - U.n.*U.vx.*U.vy;
          elseif (varid==6)  % read pressure tensor
            U.pxz= squeeze(var) - U.n.*U.vx.*U.vz;
          elseif (varid==7)  % read pressure tensor
            U.pyy= squeeze(var) - U.n.*U.vy.*U.vy;
          elseif (varid==8)  % read pressure tensor
            U.pyz= squeeze(var) - U.n.*U.vy.*U.vz;
          elseif (varid==9)  % read pressure tensor
            U.pzz= squeeze(var) - U.n.*U.vz.*U.vz;
          end
        end      
      end
      U.rest=diss(1); U.hrest=diss(2);
      fprintf('Time: %e\n',U.time)
      clear var;
      disp('Mom read');
      diary off; diary on;
      %B_spectrum_min(U,basename);       disp('B Spectrum calc')
      %V_spectrum_min(U,basename);       disp('V Spectrum calc')
      %larmor_radius(U,basename);    disp('larmor calc')
      %rateofstrain(U,basename);     disp('ROS calc')
      %wavenumbers(U,basename);      disp('Wavenumber calc')
      %PDF(U,basename);              disp('PDF calc')
      rhoBCJ(U,basename);            disp('rhoBCJ calc')
      %forcing(U,basename);          disp('Forcing calc')
    end    
  end
end

exit
