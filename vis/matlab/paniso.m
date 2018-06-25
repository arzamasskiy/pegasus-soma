
shear = 3e-4;
Lx = 1152;
Ly = 1152;

dir ='/Volumes/My Passport Studio/pegasus/fhs-slow/vtk/';
%dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/vtk/';        % directory
fname = 'fhshear';

for f=1:333

    f
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
    end

    filename = [dir,fname,'.',numlab,'.vtk'];

    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);
    timen = Grid.time;
    t(f) = timen;

    [time,var,name,status] = readvtk(Grid,filename,1);
    U.bx  = squeeze(var(1,:,:,:));
    U.by  = squeeze(var(2,:,:,:));
    U.bz  = squeeze(var(3,:,:,:));

    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d = squeeze(var);

    [time,var,name,status] = readvtk(Grid,filename,4);
    U.pxx= squeeze(var(1,:,:,:));
    U.pxz= squeeze(var(2,:,:,:));
    U.pxy= squeeze(var(3,:,:,:));
    U.pzx= squeeze(var(4,:,:,:));
    U.pzz= squeeze(var(5,:,:,:));
    U.pzy= squeeze(var(6,:,:,:));
    U.pyx= squeeze(var(7,:,:,:));
    U.pyz= squeeze(var(8,:,:,:));
    U.pyy= squeeze(var(9,:,:,:));

    bsq  = (U.bx).^2 + (U.by).^2 + (U.bz).^2;
    ptot = ( U.pxx + U.pyy + U.pzz ) / 3;
    pprl = U.bx .* ( U.pxx .* U.bx + U.pxy .* U.by + U.pxz .* U.bz ) ./ bsq ...
         + U.by .* ( U.pyx .* U.bx + U.pyy .* U.by + U.pyz .* U.bz ) ./ bsq ...
         + U.bz .* ( U.pzx .* U.bx + U.pyz .* U.by + U.pzz .* U.bz ) ./ bsq;
    pprp = 1.5 * ptot - 0.5 * pprl;
    Delta = mean(mean(1 - pprl./pprp));
    
    mun = (pprp) ./ (U.d) ./ sqrt(bsq);
    mun = log(mun);
    if (f==1)
        muo = mun;
        timeo = timen;
        nu(f) = 0;
    else
        dlnmudt = (mun-muo)/(timen - timeo);
        nu(f) = -3 * mean(mean( dlnmudt )) / Delta;
        muo = mun;
        timeo = timen;
    end
    
end