clear all;
figlab = 2;
figM = figure(figlab); hold off;

omega = 0.01;

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultTextFontSize',11);
set(0,'DefaultAxesLineWidth',1.0);


dir ='/Users/kunz/Documents/codes/pegasus/bin/mri2d-45/';        % directory
fname = 'mri';                     % filename
fn1=1;                          % first file index
fn2=67;                        % last file index
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
      U.mz = squeeze(var(2,:,:,:));
      U.my = squeeze(var(3,:,:,:));
    elseif (varid==4)    % read pressure tensor
      U.pxx= squeeze(var(1,:,:,:));
      U.pxz= squeeze(var(2,:,:,:));
      U.pxy= squeeze(var(3,:,:,:));
      U.pzx= squeeze(var(4,:,:,:));
      U.pzz= squeeze(var(5,:,:,:));
      U.pzy= squeeze(var(6,:,:,:));
      U.pyx= squeeze(var(7,:,:,:));
      U.pyz= squeeze(var(8,:,:,:));
      U.pyy= squeeze(var(9,:,:,:));
    end
    
    bxp = squeeze(mean(U.bx,2));
    %bxp = squeeze(U.bx(:,2,:));
    bxavg = mean(bxp,1);
    nfft = 2^nextpow2(Grid.nx3);
    Y = fft(bxavg,nfft)/Grid.nx3;
    Y = 2*abs(Y(1:nfft/2+1));
    yk1(i) = Y(2);
    
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
 
end

semilogy(yk1,'k');
hold on;
t = [0:1:90];
semilogy(1.15e-4*exp(0.096*t),'--k');
hold off;
axis([0 90 1e-4 1]);

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,30,30]);
set(gcf,'PaperPositionMode','auto');

%{
bxavg = mean(U.bx,1);
Grid.x2 = Grid.x2 * omega;

x2min = min(Grid.x2);
x2max = max(Grid.x2);
    
plot(Grid.x2,bxavg,'k');
hold on;
axis([x2min x2max -1 1]);
xlabel('$z \varpi_0 $');
ylabel('$\langle \delta B_x \rangle$');
set(gca,'Ticklength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.2 0.2 0.7 0.7]);
plotTickLatex2D('xtickdy',-0.01,'ytickdx',-0.01,'xlabeldy',-0.02,'ytickdy',0.005,'xtickdx',0.001);
%}

%{

figure(5); clf;

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,4,10]);
set(gcf,'PaperPositionMode','auto');

nlev = 32;
bxplot = squeeze(U.bx(:,1,:));
contourf(Grid.x1,Grid.x3,bxplot',nlev);
colormap(bluewhitered(nlev));
shading flat;

axis equal;
axis([Grid.x1min Grid.x1max Grid.x3min Grid.x3max]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('$\delta B_x$');
%xlabel('$x$');
%ylabel('$y$');
set(gca,'TickLength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.1,0.1,0.8,0.8]);

plotTickLatex2D('xtickdy',-0.005,'ytickdx',-0.01,'xlabeldy',-0.005,'ylabeldx',-0.01); 



figure(6); clf;

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,4,10]);
set(gcf,'PaperPositionMode','auto');

nlev = 32;
byplot = squeeze(U.by(:,1,:))-1/sqrt(2.);
contourf(Grid.x1,Grid.x3,byplot',nlev);
colormap(bluewhitered(nlev));
shading flat;

axis equal;
axis([Grid.x1min Grid.x1max Grid.x3min Grid.x3max]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('$\delta B_y$');
%xlabel('$x$');
%ylabel('$y$');
set(gca,'TickLength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.1,0.1,0.8,0.8]);

plotTickLatex2D('xtickdy',-0.005,'ytickdx',-0.01,'xlabeldy',-0.005,'ylabeldx',-0.01); 

%}


