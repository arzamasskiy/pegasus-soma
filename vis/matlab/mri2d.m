clear all;
figlab = 1;
figM = figure(figlab); hold off;

omega = 0.01;

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
set(0,'DefaultAxesLineWidth',1.0);


dir ='/Users/kunz/Documents/codes/pegasus/bin/betaz10_ang45/';        % directory
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
      U.bz = squeeze(var(2,:,:,:));
      U.by = squeeze(var(3,:,:,:));
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
    
    %{
    bxavg = mean(U.bx,1);
    byavg = mean(U.by,1)-1/sqrt(2);
    bzavg = mean(U.bz,1)-1/sqrt(2);
    
    
    nfft = 2^nextpow2(Grid.nx2);
    Y = fft(bxavg,nfft)/Grid.nx2;
    Y = (2*abs(Y(1:nfft/2+1))).^2;
    yk1(i+1) = Y(2);
    yk2(i+1) = Y(3);
    yk3(i+1) = Y(4);
    yk4(i+1) = Y(5);
    yk(i+1) = sum(Y);
    t(i+1) = Grid.time;
    %}
    
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

%{
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,10,6]);
set(gcf,'PaperPositionMode','auto')

orb = t*omega/(2*pi);
semilogy(orb,yk,'k');
hold on;
semilogy(orb,yk1,'b');
semilogy(orb,yk2,'r');
%semilogy(orb,yk3,'g');
%semilogy(orb,yk4,'c');
leg = legend('$~\sum_k$','$~k=1$','$~k=2$','Location','Southeast');
legend(leg,'boxoff');
tt = 0:0.1:2;
%
% angle = 45, beta = 10
%
semilogy(tt,9e-9*exp(2*0.96*2*pi*tt),'--b');
semilogy(tt,2.8e-7*exp(2*0.492*2*pi*tt),'--r');
%semilogy(tt,2.3e-7*exp(-2*0.2886*2*pi*tt),'--g');
%
%
%
hold off;
axis([0 1.5 1e-9 1]);
xlabel('time in orbits');
ylabel('$| \delta B_x |^2$');
set(gca,'YTick',[1e-8 1e-6 1e-4 1e-2 1e0]);
set(gca,'TickLength',[0.02 0.02]);
set(gca,'XMinorTick','on')
set(gca,'Units','normalized','Position',[0.17,0.15,0.77,0.82]);

plotTickLatex2D('xtickdy',-0.01,'xtickdx',0.005,'ytickdx',-0.01,'ylabeldx',-0.01,'xlabeldy',-0.01);
%}

%{

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,10,6]);
set(gcf,'PaperPositionMode','auto');

Grid.x2 = Grid.x2 * omega;
x2min = min(Grid.x2);
x2max = max(Grid.x2);
    
plot(Grid.x2,bxavg,'r');
hold on;
plot(Grid.x2,byavg,'b');
plot(Grid.x2,bzavg,'k');
axis([x2min x2max -1 1]);
xlabel('$z \varpi_0 $');
set(gca,'Ticklength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.16,0.15,0.76,0.82]);
plotTickLatex2D('xtickdy',-0.01,'ytickdx',-0.005,'xlabeldy',-0.03,'ytickdy',0.005,'xtickdx',0.005);


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

%}


for j=1:Grid.nx2
	xx(:,j) = Grid.x1(:);
end
B0 = 1/sqrt(2);
dbz = U.bz-B0;
dbx = U.bx;

Lx = 1.430597747536336e+02;
Lz = 5.722390990145344e+02;
kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
kx = kx1 * (2*pi/Grid.dx1);
kz1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
kz = kz1 * (2*pi/Grid.dx2);
[KX,KZ] = meshgrid(kx,kz);
K2 = KX.^2 + KZ.^2;
K2(1,1) = 1.0;
ay = ifft2( 1i./K2' .* ( KZ' .* fft2(dbx) - KX' .* fft2(dbz) ) );
ay = real(ay) + xx*B0;


ax1 = axes;
imagesc(Grid.x1,Grid.x2,U.bx');
axis image; axis xy;
ax2 = axes('Position',get(ax1,'Position'),...
           'Box','off','Color','none',...
           'Ylim',get(ax1,'YLim'),'Xlim',get(ax1,'XLim'));
hold on;
contour(ax2,Grid.x1,Grid.x2,ay',24,'-k');
axis image; axis xy;
hold off;
