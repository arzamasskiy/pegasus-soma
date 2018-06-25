clear all;

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
set(0,'DefaultAxesLineWidth',1.0);


dir ='/Users/kunz/';        % directory
fname = 'shock';                     % filename
fn1=60;                          % first file index
fn2=60;                        % last file index

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

    for varid=1:4
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
    end
    
    bxavg = fliplr(mean(U.bx,2));
    byavg = fliplr(mean(U.by,2));
    bzavg = fliplr(mean(U.bz,2));
    davg  = fliplr(mean(U.d ,2));
    
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

for i=1:Grid.nx1
    for j=1:Grid.nx2
        dBy(i+Grid.nx1,j) = U.by(i,j);
        dBy(i,j) = U.by(Grid.nx1-i+1,j);
        dBx(i+Grid.nx1,j) = U.bx(i,j)-1.0;
        dBx(i,j) = U.bx(Grid.nx1-i+1,j)-1.0;
    end
end

for i=1:2*Grid.nx1
    yy(i,:) = Grid.x2(:);
end

kx1 = mod( 1/2 + (0:(2*Grid.nx1-1))/(2*Grid.nx1) , 1 ) - 1/2;
kx = kx1 * (2*pi/Grid.dx1);
ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
ky = ky1 * (2*pi/Grid.dx2);
[KX,KY] = meshgrid(kx,ky);
K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
az2 = ifft2( 1i./K2' .* ( KY' .* fft2(dBx) - KX' .* fft2(dBy) ) );
az2 = real(az2) + yy;
az(1:Grid.nx1,:) = az2(Grid.nx1+1:2*Grid.nx1,:);
az = fliplr(az);
clear az2; clear dBy; clear dBx;


fig(1);clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,16.5,12]);
set(gcf,'PaperPositionMode','auto');

%den = U.d - min(min(U.d));
%den = den/max(max(den));
%den = den - 0.5;
den = fliplr(U.d);
bz = fliplr(U.bz);

x = fliplr(Grid.x1);

subplot(3,1,1);
imagesc(x,Grid.x2,den',[0.8 5]);
colormap bluewhitered;
axis xy; shading flat;
%freezeColors;
ylabel('$y$');
set(gca,'XTickLabel',[]);
set(gca,'YTick',[20 40 60 80]);
set(gca,'TickLength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.08,0.71,0.88,0.28]);
plotTickLatex2D('ytickdx',0.005,'ylabeldx',-0.01,'ytickdy',0.007);

subplot(3,1,2);
imagesc(x,Grid.x2,bz',[-2 2]);
colormap hawley;
axis xy; shading flat;
ylabel('$y$');
set(gca,'XTickLabel',[]);
set(gca,'YTick',[20 40 60 80]);
set(gca,'TickLength',[0.02 0.02]);
set(gca,'Units','normalized','Position',[0.08,0.40,0.88,0.28]);
plotTickLatex2D('ytickdx',0.005,'ylabeldx',-0.01,'ytickdy',0.007);


subplot(3,1,3);
contour(x,Grid.x2,az',20,'k');
xlabel('$x$'); ylabel('$y$');
set(gca,'TickLength',[0.02 0.02]);
set(gca,'YTick',[20 40 60 80]);
set(gca,'XTick',[0 500 1000 1500 2000]);
set(gca,'Units','normalized','Position',[0.08,0.09,0.88,0.28]);
plotTickLatex2D('ytickdx',0.005,'ylabeldx',-0.01,'ytickdy',0.007,'xlabeldy',0.005);


fig(3);clf;
plot(x,davg);
axis([0 2048 0 5]);

fig(4);clf;
plot(x,bzavg);
axis([0 2048 -1.5 1.5]);


