clear all;
figlab = 1;
figM = figure(figlab); hold off;

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultTextFontSize',11);
set(0,'DefaultAxesLineWidth',1.0);

dir ='/Users/kunz/Documents/codes/pegasus/bin/mri2d-znf/hyperres/secondtry/';        % directory
fname = 'mri2d';                     % filename
fn1=90;                          % first file index
fn2=90;                        % last file index
varid=1;                        % variable id: 1 = b, 2 = d, 3 = m, 4 = p

for i=fn1:fn2

    % format file number
    if (i<10)
      numlab = ['000',num2str(i)];
    elseif (i<100)
      numlab = ['00',num2str(i)];
    elseif (i<1000)
      numlab = ['0',num2str(i)];
    else
      numlab = num2str(i);
    end
 
    % declare file name
    filename = [dir,fname,'.',numlab,'.vtk'];

    % open file and initialize grid
    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);

    [time,var,name,status] = readvtk(Grid,filename,varid);
    if (varid==1)         % read magnetic field
      U.bx = squeeze(var(1,:,:,:));
      U.by = squeeze(var(2,:,:,:));
      U.bz = squeeze(var(3,:,:,:));
    elseif (varid==2)    % read density
      U.d = squeeze(var);
    elseif (varid==3)    % read momentum density
      U.mx = squeeze(var(1,:,:,:));
      U.my = squeeze(var(2,:,:,:));
      U.mz = squeeze(var(3,:,:,:));
    elseif (varid==4)    % read pressure tensor
      U.pxx= squeeze(var(1,:,:,:));
      U.pxy= squeeze(var(2,:,:,:));
      U.pxz= squeeze(var(3,:,:,:));
      U.pyx= squeeze(var(4,:,:,:));
      U.pyy= squeeze(var(5,:,:,:));
      U.pyz= squeeze(var(6,:,:,:));
      U.pzx= squeeze(var(7,:,:,:));
      U.pzy= squeeze(var(8,:,:,:));
      U.pzz= squeeze(var(9,:,:,:));
    end
    
    % pressure anisotropy
    %{
      bsq  = (U.bx).^2 + (U.by).^2 + (U.bz).^2;
      ptot = ( U.pxx + U.pyy + U.pzz ) / 3;
      pprl = U.bx .* ( U.pxx .* U.bx + U.pxy .* U.by + U.pxz .* U.bz ) ./ bsq ...
           + U.by .* ( U.pyx .* U.bx + U.pyy .* U.by + U.pyz .* U.bz ) ./ bsq ...
           + U.bz .* ( U.pzx .* U.bx + U.pyz .* U.by + U.pzz .* U.bz ) ./ bsq;
      pprp = 1.5 * ptot - 0.5 * pprl;
      Delta= 1.0 - pprl./pprp;
    %}

    
   % plot(Grid.x1,U.bz,'b');
    %hold on;
    %plot(Grid.x1,U.by,'r');
    %plot(Grid.x1,Delta,'k');
    %hold off;
    %axis([-150 150 -0.000015 0.000015]);
    %getframe(figM);
    %}
 
end

%{
figure(5); clf;

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,10,10]);
set(gcf,'PaperPositionMode','auto');

nlev = 32;
contourf(Grid.x1,Grid.x2,U.by',nlev);
colormap(bluewhitered(nlev));
shading flat;

axis equal;
axis([Grid.x1min Grid.x1max Grid.x1min Grid.x2max]);
%set(gca,'XTick',[0 1 2 3 4]);
%set(gca,'YTick',[0.5 1 1.5 2]);
xlabel('$x$');
ylabel('$y$');
set(gca,'TickLength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.2,0.2,0.7,0.7]);

plotTickLatex2D('xtickdy',-0.005,'ytickdx',-0.01,'xlabeldy',-0.005,'ylabeldx',-0.01); 

%}