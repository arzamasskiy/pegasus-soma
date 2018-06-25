set(0,'DefaultLineLineWidth',2.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',40);
set(0,'DefaultTextFontSize',40);
set(0,'DefaultAxesLineWidth',2.0);

 


fn1=0;
fn2=628;
%track1 = 900;
%track2 = 99;

shear = 3e-4; %3e-4; %1e-4;
Lx = 1152;
Ly = 1152;

dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/';        % directory
fname = 'vtk/mrshear';
%pname = 'track.lis/mrshear';

for ff=fn1:fn2
    
    i = ff-fn1+1;

    
    % format file number
    if (ff<10)
        numlab = ['000',num2str(ff)];   
    elseif (ff<100)
        numlab = ['00',num2str(ff)];
    elseif (ff<1000)
        numlab = ['0',num2str(ff)];
    else
        numlab = num2str(ff);
    end


    % declare file name
    filename = [dir,fname,'.',numlab,'.vtk'];

    % open file and initialize grid
    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);

    [time,var,name,status] = readvtk(Grid,filename,1);
    U.bx  = squeeze(var(1,:,:,:));
    U.by  = squeeze(var(2,:,:,:));
    U.bz  = squeeze(var(3,:,:,:));

    Bsq = U.bx.^2 + U.by.^2 + U.bz.^2;
    
    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d = squeeze(var);
    
    strain = ((Bsq).^(1.5)) ./ ((U.d).^(2));
    
    avg(i) = sum(sum(strain))/Grid.nx1/Grid.nx2;
    
end

dt = 10;
logavg = log(avg);
strain = diff(logavg)/dt;

dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/';        % directory
fname = 'mrshear';
filename = [dir,fname,'.hst'];

fid = fopen(filename);
C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);

time = C{1};
Delt = C{22};

for i=1:fn2-fn1+1
    t(i) = time((i-1)*dt+1);
    paniso(i) = Delt((i-1)*dt+1);
end

Del = paniso(1:628); time = t(1:628);
dlnDdt = diff(paniso)./Del/dt;

nu = 3*strain./Del - dlnDdt;

plot(t*shear,nu);
