clear all;

i=10;

%1 = contour plots for FIREHOSE
%2 = contour plots for MIRROR
%3 = box average evolution for FIREHOSE
%4 = box average evolution for MIRROR
%5 = spectrum for MIRROR
%6 = spectrum for FIREHOSE
%7 = distribution function for MIRROR
%8 = distribution function for FIREHOSE
%9 = particle tracks for MIRROR
%10 = particle tracks for FIREHOSE
%11 = nueff plot for FIREHOSE & MIRROR
%55 = saturated spectra for MIRROR & FIREHOSE

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
set(0,'DefaultAxesLineWidth',1.0);

format long;

if i==1
    
  figure(1);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,5.5]);
  set(gcf,'PaperPositionMode','auto');
 
  Lx = 1152;
  Ly = 1152;
  
  shear = 3e-4;
  dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/vtk/';        % directory
  fname = 'fhshear';
  
  % contour limits
  clims = [-0.8 0.8];
  % field lines (isocontours of vector potential)
  alin = linspace(-1.1273e+03,1.1162e+03,30);

%{
  subplot(1,2,1);
  f=22;
  
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

  [time,var,name,status] = readvtk(Grid,filename,1);
  U.bx  = squeeze(var(1,:,:,:));
  U.by  = squeeze(var(2,:,:,:));
  U.bz  = squeeze(var(3,:,:,:));
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  set(gca,'XTickLabel',[' ']); ylabel('$y$');title('$St = 0.066$');
  set(gca,'Units','normalized','Position',[0.11,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);

  subplot(1,2,2);
  f=50;
  
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

  [time,var,name,status] = readvtk(Grid,filename,1);
  U.bx  = squeeze(var(1,:,:,:));
  U.by  = squeeze(var(2,:,:,:));
  U.bz  = squeeze(var(3,:,:,:));
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  set(gca,'YTickLabel',[' '],'XTickLabel',[' ']);title('$St = 0.15$');
  set(gca,'Units','normalized','Position',[0.504,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);

  cb = colorbar;
  set(cb,'Position',[0.894,0.14,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-0.8; ymax=0.8;
  set(cb,'YTick',ymin:0.2:ymax);
  yt_label={'$-0.8$','$-0.6$','$-0.4$','$-0.2$','$0.0$','$0.2$','$0.4$','$0.6$','$0.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
  figure(10);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,5.5]);
  set(gcf,'PaperPositionMode','auto');
 
  
  subplot(1,2,1);
  f=333;
  
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

  [time,var,name,status] = readvtk(Grid,filename,1);
  U.bx  = squeeze(var(1,:,:,:));
  U.by  = squeeze(var(2,:,:,:));
  U.bz  = squeeze(var(3,:,:,:));
  
  [time,var,name,status] = readvtk(Grid,filename,2);
  U.d   = squeeze(var);
  U.d   = U.d - mean(mean(U.d));
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x$'); ylabel('$y$'); title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.11,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);

  subplot(1,2,2);
  
  imagesc(Grid.x1,Grid.x2,220*U.d',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  set(gca,'TickLength',[0.02 0.02]);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x$'); set(gca,'YTickLabel',[' ']);title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.504 0.08 0.365 0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);
  annotation('textbox',[0.76 0.159 0.1 0.09],'FontSize',10, ...
    'Interpreter','latex','String','2$\hspace{-0.4mm}$2$\hspace{-0.2mm}$0$\delta\!n$','Margin',3, ...
    'BackgroundColor','w','EdgeColor','k','VerticalAlignment','baseline');

  cb = colorbar;
  set(cb,'Position',[0.894,0.14,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-0.8; ymax=0.8;
  set(cb,'YTick',ymin:0.2:ymax);
  yt_label={'$-0.8$','$-0.6$','$-0.4$','$-0.2$','$0.0$','$0.2$','$0.4$','$0.6$','$0.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
  %}
  
  
  %{
  subplot(1,2,1);
  f=167;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  %dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);

  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x$'); ylabel('$y$');title('$St = 0.5$');
  set(gca,'Units','normalized','Position',[0.11,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);

  subplot(1,2,2);
  f=333;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
  
  [time,var,name,status] = readvtk(Grid,filename,2);
  U.d   = squeeze(var);
  U.d   = U.d - mean(mean(U.d));
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  %dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x$'); set(gca,'YTickLabel',[' ']);title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.504,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);
%  annotation('textbox',[0.76 0.159 0.1 0.09],'FontSize',10, ...
%    'Interpreter','latex','String','2$\hspace{-0.4mm}$7$\hspace{-0.2mm}$0$\delta\!n$','Margin',3, ...
%    'BackgroundColor','w','EdgeColor','k','VerticalAlignment','baseline');

  cb = colorbar;
  set(cb,'Position',[0.894,0.14,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-0.8; ymax=0.8;
  set(cb,'YTick',ymin:0.2:ymax);
  yt_label={'$-0.8$','$-0.6$','$-0.4$','$-0.2$','$0.0$','$0.2$','$0.4$','$0.6$','$0.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
    %}
  
  
  
  subplot(1,2,1);
  f=22;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x/d_{\!\rm i0}$'); ylabel('$y/d_{\!\rm i0}$');title('$St = 0.066$');
  set(gca,'Units','normalized','Position',[0.11,0.095,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035,'ylabeldx',0.011);

  subplot(1,2,2);
  
  f=333;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
  
  [time,var,name,status] = readvtk(Grid,filename,2);
  U.d   = squeeze(var);
  U.d   = U.d - mean(mean(U.d));
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = 3/sqrt(13);
  bx0 = 2/sqrt(13);
  by  = by0-bx0*shear*time;
  dBz = U.bz;
  dBy = U.by-by;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by*xx + az;
  
  imagesc(Grid.x1,Grid.x2,dBz',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1.1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x/d_{\!\rm i0}$'); set(gca,'YTickLabel',[' ']); title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.504 0.095 0.365 0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035);
  
  cb = colorbar;
  set(cb,'Position',[0.894,0.155,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-0.8; ymax=0.8;
  set(cb,'YTick',ymin:0.2:ymax);
  yt_label={'$-0.8$','$-0.6$','$-0.4$','$-0.2$','$0.0$','$0.2$','$0.4$','$0.6$','$0.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
  
  

elseif i==2
%
% MIRROR CONTOUR PLOTS
%

  figure(2);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,5.5]);
  set(gcf,'PaperPositionMode','auto');
 
  Lx = 1152;
  Ly = 1152;
  
  shear = 3e-4;
  dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/vtk/';
  %shear = 1e-4;
  %dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-super-slow/vtk/';
  fname = 'mrshear';
  
  clims = [-1.8 1.8];
  alin = linspace(-1.49e+03,1.56e+03,30);
  
  subplot(1,2,1);
  f=50;
      
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = (-1/sqrt(5)-2/sqrt(5)*shear*time);
  bx0 = 2/sqrt(5);
  dBz = U.bz;
  dBy = U.by-by0;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by0*xx + az;
  
  dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);

  imagesc(Grid.x1,Grid.x2,dBprl',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1);
  hold off;
  set(gca,'TickLength',[0.02 0.02],'FontSize',10);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  set(gca,'XTickLabel',[' ']); ylabel('$y/d_{\!\rm i0}$');title('$St = 0.15$');
  set(gca,'Units','normalized','Position',[0.11,0.095,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035,'ylabeldx',0.011);

  subplot(1,2,2);
  f=166;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = (-1/sqrt(5)-2/sqrt(5)*shear*time);
  bx0 = 2/sqrt(5);
  dBz = U.bz;
  dBy = U.by-by0;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by0*xx + az;
  
  dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);
  
  
  imagesc(Grid.x1,Grid.x2,dBprl',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  set(gca,'YTickLabel',[' '],'XTickLabel',[' ']); title('$St = 0.5$');
  set(gca,'Units','normalized','Position',[0.504,0.095,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035);


  cb = colorbar;
  set(cb,'Position',[0.894,0.155,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-1.8; ymax=1.8;
  set(cb,'YTick',ymin:0.4:ymax);
  yt_label={'$-1.8$','$-1.4$','$-1.0$','$-0.6$',...
      '$-0.2$','$0.2$','$0.6$','$1.0$',...
      '$1.4$','$1.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end

  
  figure(20);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,5.5]);
  set(gcf,'PaperPositionMode','auto');

  subplot(1,2,1);
  f=333;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
  end

  % declare file name
  filename = [dir,fname,'.',numlab,'.vtk'];

  % open file and initialize grid
  % B, d, M, P
  ary = [3 1 3 9];
  [Grid,status] = init_grid(filename,ary);

  [time,var,name,status] = readvtk(Grid,filename,1);
  U.bx  = squeeze(var(1,:,:,:));
  U.by  = squeeze(var(2,:,:,:));
  U.bz  = squeeze(var(3,:,:,:));
  
  [time,var,name,status] = readvtk(Grid,filename,2);
  U.d = squeeze(var);
  U.d = U.d - mean(mean(U.d));
    
  %{
  [time,var,name,status] = readvtk(Grid,filename,3);
  U.vx = squeeze(var(1,:,:,:))./U.d;
  U.vy = squeeze(var(2,:,:,:))./U.d;
%}
  
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = (-1/sqrt(5)-2/sqrt(5)*shear*time);
  bx0 = 2/sqrt(5);
  dBz = U.bz;
  dBy = U.by-by0;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by0*xx + az;
  
  dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);

  imagesc(Grid.x1,Grid.x2,dBprl',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1);
  hold off;
  set(gca,'TickLength',[0.02 0.02],'FontSize',10);
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x/d_{\!\rm i0}$'); ylabel('$y/d_{\!\rm i0}$');title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.11,0.095,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035,'ylabeldx',0.011);

  subplot(1,2,2);
  
  imagesc(Grid.x1,Grid.x2,130*U.d',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  %{
  hold on;
  imax=36;
  skip=1152/imax;
for i=1:imax
    xvx(i) = Grid.x1(i*skip);
    yvy(i) = Grid.x2(i*skip);
    for j=1:imax
        vx(i,j) = U.vx(i*skip,j*skip);
        vy(i,j) = U.vy(i*skip,j*skip);
    end
end
  quiver(xvx,yvy,vx',vy',2,'k','Linewidth',0.4);
  hold off;
  %}
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x/d_{\!\rm i0}$'); set(gca,'YTickLabel',[' ']);title('$St = 1.0$');
  set(gca,'Units','normalized','Position',[0.504,0.095,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.035);
  annotation('textbox',[0.75 0.174 0.11 0.09],'FontSize',10, ...
    'Interpreter','latex','String','1$\hspace{-0.4mm}$3$\hspace{-0.2mm}$0$\delta\!n_{\!{\rm i}}$','Margin',3, ...
    'BackgroundColor','w','EdgeColor','k','VerticalAlignment','baseline');

  cb = colorbar;
  set(cb,'Position',[0.894,0.155,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-1.8; ymax=1.8;
  set(cb,'YTick',ymin:0.4:ymax);
  yt_label={'$-1.8$','$-1.4$','$-1.0$','$-0.6$',...
      '$-0.2$','$0.2$','$0.6$','$1.0$',...
      '$1.4$','$1.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  

  
  
  
  %{
  f=490;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = (-1/sqrt(5)-2/sqrt(5)*shear*time);
  bx0 = 2/sqrt(5);
  dBz = U.bz;
  dBy = U.by-by0;
  dBx = U.bx-bx0;

  dBxp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBx);
  dByp=remap2d(Grid.x1,Grid.x2,shear,time,-1,dBy);
  kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
  kx = kx1 * (2*pi/Grid.dx1);
  ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
  ky = ky1 * (2*pi/Grid.dx2);
  tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
  [KX,KY] = meshgrid(kx,ky);
  KX = KX + KY*shear*tremap;
  K2 = KX.^2 + KY.^2; K2(1,1) = 1.0;
  azp = ifft2( 1i./K2' .* ( - KY' .* fft2(dBxp) + KX' .* fft2(dByp) ) );
  azp = real(azp);
  az = remap2d(Grid.x1,Grid.x2,shear,time,1,azp);
  az = bx0*yy - by0*xx + az;
  
  dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);
  
  
  imagesc(Grid.x1,Grid.x2,dBprl',clims);
  shading flat; axis xy; axis image;
  colormap hawley;
  hold on;
  contour(Grid.x1,Grid.x2,az',alin,'k','LineWidth',1);
  hold off;
  set(gca,'TickLength',[0.02 0.02]); 
  set(gca,'YTick',[-500:250:500]); set(gca,'XTick',[-500:250:500]);
  xlabel('$x$'); set(gca,'YTickLabel',[' ']);title('$St = 1.5$');
  set(gca,'Units','normalized','Position',[0.504,0.08,0.365,0.88]);
  plotTickLatex2D('ytickdy',0.005,'xtickdy',-0.02,'xlabeldy',-0.04);


  cb = colorbar;
  set(cb,'Position',[0.894,0.14,0.025,0.76],'Ticklength',[0.015 0.015]);
  ymin=-1.8; ymax=1.8;
  set(cb,'YTick',ymin:0.4:ymax);
  yt_label={'$-1.8$','$-1.4$','$-1.0$','$-0.6$',...
      '$-0.2$','$0.2$','$0.6$','$1.0$',...
      '$1.4$','$1.8$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  %}
  
  
  %{
  f=166;
  
  % format file number
  if (f<10)
    numlab = ['000',num2str(f)];   
  elseif (f<100)
    numlab = ['00',num2str(f)];
  elseif (f<1000)
    numlab = ['0',num2str(f)];
  else
    numlab = num2str(f);
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
    
  [time,var,name,status] = readvtk(Grid,filename,2);
  U.d = squeeze(var);
    
  for j=1:Grid.nx1
    yy(j,:) = Grid.x2(:);
  end
  for i=1:Grid.nx2
    xx(:,i) = Grid.x1(:);
  end
  by0 = (-1/sqrt(5)-2/sqrt(5)*shear*time);
  bx0 = 2/sqrt(5);
  dBz = U.bz;
  dBy = U.by-by0;
  dBx = U.bx-bx0;
  dBprl = (by0*dBy+bx0*dBx)/sqrt(by0*by0+bx0*bx0);

  

  fig(99);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,9]);
  set(gcf,'PaperPositionMode','auto');
  
  gax = axes('Position',[0.01 0.58 0.35 0.35]);    
  clims = [0 3];%[-1 1];
  imagesc(Grid.x1,Grid.x2,100*(U.d'-1),clims);
  shading flat; axis xy; axis image;
  colormap bluewhitered;
  set(gca,'TickLength',[0.02 0.02],'FontSize',11);
  set(gca,'YTick',[-400:200:400]); set(gca,'XTick',[-400:200:400]);
  xlabel('$x$'); ylabel('$y$');
  set(gca,'Units','normalized','Position',[0.065,0.1,0.82,0.88]);

  plotTickLatex2D('ytickdy',0.006,'ytickdx',-0.005);

  cb = colorbar;
  set(cb,'Position',[0.86,0.1,0.04,0.88],'Ticklength',[0.015 0.015]);
  ymin=-0.13; ymax=3.13;
  set(cb,'YTick',ymin:0.5:ymax);
  yt_label={'$-1.0$','$-0.8$','$-0.6$','$-0.4$','$-0.2$','$0$','$0.2$','$0.4$','$0.6$','$0.8$','$1.0$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
  fig(100);clf;
  set(gcf,'Color',[1 1 1]);
  set(gcf,'Units','centimeters');
  set(gcf,'Position',[1,1,11.5,9]);
  set(gcf,'PaperPositionMode','auto');
  
  gax = axes('Position',[0.01 0.58 0.35 0.35]);    
  clims = [-1.6 1.6];
  imagesc(Grid.x1,Grid.x2,dBprl',clims);
  shading flat; axis xy; axis image;
  colormap bluewhitered;
  set(gca,'TickLength',[0.02 0.02],'FontSize',11);
  set(gca,'YTick',[-400:200:400]); set(gca,'XTick',[-400:200:400]);
  xlabel('$x$'); ylabel('$y$');
  set(gca,'Units','normalized','Position',[0.065,0.1,0.82,0.88]);

  plotTickLatex2D('ytickdy',0.006,'ytickdx',-0.005);

  cb = colorbar;
  set(cb,'Position',[0.86,0.1,0.04,0.88],'Ticklength',[0.015 0.015]);
  ymin=-1.6; ymax=1.6;
  set(cb,'YTick',ymin:0.4:ymax);
  yt_label={'$-1.6$','$-1.2$','$-0.8$','$-0.4$','$0$','$0.4$','$0.8$','$1.2$','$1.6$'};
  h=gcf; c=get(h,'children');
  set(h,'CurrentAxes',cb);
  set(cb,'YTickLabel',[]);
  yt_pos=linspace(ymin,ymax,(size(yt_label,2)));
  for k=1:size(yt_label,2)
    t(k)=text(3,yt_pos(k),yt_label{k},'HorizontalAlignment','left');
    set(t(k),'FontSize',10);
  end
  
  %}
  

elseif i==3
%
% FIREHOSE
% box average evolution of magnetic energy; Reynolds/Maxwell/pressure 
% stresses; instability parameter
%
    shear = 3e-3;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-fast/';        % directory
    fname = 'fhshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,repmat('%f',[1,24]),'CommentStyle','#');
    fclose(fid);
    
    
    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam1 = C{21};
    Del1 = C{22};
    
    time1 = time*shear;

    ME1 = MEx+MEy+MEz;
    B1sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time1).^2);
    by = (pitch-time1)./sqrt(1.0+(pitch-time1).^2);
    dBprlsq1 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq1 = 2*ME1-B1sq-dBprlsq1;
    

    beta = 100./ME1;   
    for i=1:length(time1)-1
        dlnBdt1(i) = 0.5*(log(ME1(i+1))-log(ME1(i)));
    end
    dlnBdt1(i+1) = dlnBdt1(i);
    num1 = -1.5.*beta(:).*dlnBdt1(:); %shear*beta.*bx.*by;
    
        
    shear = 1e-3;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-medium/';        % directory
    fname = 'fhshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam2 = C{21};
    Del2 = C{22};
    
    time2 = time*shear;

    ME2 = MEx+MEy+MEz;
    B2sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time2).^2);
    by = (pitch-time2)./sqrt(1.0+(pitch-time2).^2);
    dBprlsq2 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq2 = 2*ME2-B2sq-dBprlsq2;
       

    beta = 100./ME2;   
    for i=1:length(time2)-1
        dlnBdt2(i) = 0.5*(log(ME2(i+1))-log(ME2(i)));
    end
    dlnBdt2(i+1) = dlnBdt2(i);
    num2 = -1.5.*beta(:).*dlnBdt2(:); %shear*beta.*bx.*by;
    
    
    shear = 3e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/';        % directory
    fname = 'fhshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam3 = C{21};
    Del3 = C{22};
    
    time3 = time*shear;

    ME3 = MEx+MEy+MEz;
    B3sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time3).^2);
    by = (pitch-time3)./sqrt(1.0+(pitch-time3).^2);
    dBprlsq3 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq3 = 2*ME3-B3sq-dBprlsq3;

    beta = 100./ME3;   
    for i=1:length(time3)-1
        dlnBdt3(i) = 0.5*(log(ME3(i+1))-log(ME3(i)));
    end
    dlnBdt3(i+1) = dlnBdt3(i);
    num3 = -1.5.*beta(:).*dlnBdt3(:); %shear*beta.*bx.*by;
    
    
    shear = 1e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-super-slow/';        % directory
    fname = 'fhshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam4 = C{21};
    Del4 = C{22};
    
    time4 = time*shear;

    ME4 = MEx+MEy+MEz;
    B4sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time4).^2);
    by = (pitch-time4)./sqrt(1.0+(pitch-time4).^2);
    dBprlsq4 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya); 
    dBprpsq4 = 2*ME4-B4sq-dBprlsq4;
    
    beta = 100./ME4;   
    for i=1:length(time4)-1
        dlnBdt4(i) = 0.5*(log(ME4(i+1))-log(ME4(i)));
    end
    dlnBdt4(i+1) = dlnBdt4(i);
    num4 = -1.5*beta(:).*dlnBdt4(:); %shear*beta.*bx.*by;
    
    %{
    shear = 3e-5;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-super-duper-slow/';        % directory
    fname = 'fhshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam5 = C{21};
    Del5 = C{22};
    
    time5 = time*shear;

    ME5 = MEx+MEy+MEz;
    B5sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time5).^2);
    by = (pitch-time5)./sqrt(1.0+(pitch-time5).^2);
    dBprlsq5 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq5 = 2*ME5-B5sq-dBprlsq5;

    beta = 200./B5sq;   
    num5 = shear*200;
 %}
    
    
    
    figure(3);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    
    %
    % if you want one panel with dBprp
    %  
    set(gcf,'Position',[1,1,11.5,5.5]);
    
    h(1) = loglog(time1,dBprpsq1,'Color',[232/255 0 0]);
    hold on;
    h(2) = loglog(time2,dBprpsq2,'Color',[34/255 139/255 34/255]);
    h(3) = loglog(time3,dBprpsq3,'b');
    h(4) = loglog(time4,dBprpsq4,'k');
    t = linspace(3e-2,4.2e-1,2);
    loglog(t,t*1.2,'--k','LineWidth',.4);
    hold off;
    axis([2e-2 3 1e-7 1e0]);
    xlabel('$S t$'); ylabel('$\langle |\delta\!B_{\!\perp}|^{2} \rangle / B^2_0$');
    set(gca,'XTick',[1e-1 1e0]);
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.83]);
    plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.007,'ylabeldx',-0.025,'ytickdy',0.004);
    text(7e-2,2.6e-1,'$\propto\! t$','FontSize',10,'Interpreter','latex');
    text(0.025,0.15,'(a)','Fontsize',10);
    text(2.5e-2,1e-4,'$S = 10^{-4}$','Fontsize',8,'Rotation',74);
    text(4.85e-2,1.8e-4,'$3 \times 10^{-4}$','Fontsize',8,'Rotation',73.5);
    text(1.035e-1,6e-4,'$S = 10^{-3}$','Fontsize',8,'Rotation',71.9);
    text(2.08e-1,8e-4,'$3 \times 10^{-3}$','Fontsize',8,'Rotation',71);
    
    %{
    hax = axes('Position',[0.65 0.34 0.28 0.37]);
    S = [1e-4,3e-4,1e-3,3e-3];
    BL = linspace(-4.2,-2.35,2); BD = [0.07204,.125,.23,.3366]; 
    plot(BL,0.5*BL+0.86,'--k','LineWidth',0.5);
    hold on;
    plot(log10(S(1)),log10(BD(1)),'+k');
    plot(log10(S(2)),log10(BD(2)),'+b');
    plot(log10(S(3)),log10(BD(3)),'+','Color',[34/255 139/255 34/255]);
    plot(log10(S(4)),log10(BD(4)),'+','Color',[232/255 0 0]);
    
    axis([-4.5 -2 -1.4 -0.2]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10);
    set(gca,'XMinorTick','on','YMinorTick','on','YTick',[-1 -0.5]);
    xlabel('${\rm lg}~S$'); title('${\rm lg}~\langle | \delta\!B^2_\perp | \rangle_{\rm sat}$');
    plotTickLatex2D('xlabeldy',-0.01,'xtickdy',-0.01,'ytickdx',-0.005,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(-4,-0.52,'$\propto\!S^{1/2}$','FontSize',10,'Interpreter','latex');    
    hold off;
    %}
    hax = axes('Position',[0.62 0.34 0.3 0.37]);
    S = [1e-4,3e-4,1e-3,3e-3];
    BL = linspace(-4.15,-2.35,2); BD = [0.07204,.125,.23,.3366]; 
    loglog(10.^BL,10.^(0.5*BL+0.86),'--k','LineWidth',0.5);
    hold on;
    loglog((S(1)),(BD(1)),'+k');
    loglog((S(2)),(BD(2)),'+b');
    loglog((S(3)),(BD(3)),'+','Color',[34/255 139/255 34/255]);
    loglog((S(4)),(BD(4)),'+','Color',[232/255 0 0]);
    axis([3e-5 1e-2 0.04 0.6]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10,'XTick',[1e-4 1e-3 1e-2]);
    set(gca,'XMinorTick','on','YMinorTick','on','YTick',10.^[-1]);
    xlabel('$S$'); title('$\langle | \delta\!B_\perp |^2 \rangle_{\rm sat}~/ B^2_0$');
    plotTickLatex2D('xlabeldy',-0.01,'xtickdy',-0.01,'ytickdx',-0.002,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(1.2e-4,10^-0.52,'$\propto\!S^{1/2}$','FontSize',10,'Interpreter','latex');    
    hold off;
    
    %
    % if you want two panels with dBprp and dBprl
    %
    %{
    set(gcf,'Position',[1,1,11.5,10.5]);
    
    subplot(2,1,1);
    
    
    h(1) = loglog(time1,dBprpsq1,'Color',[232/255 0 0]);
    hold on;
    h(2) = loglog(time2,dBprpsq2,'Color',[34/255 139/255 34/255]);
    h(3) = loglog(time3,dBprpsq3,'b');
    h(4) = loglog(time4,dBprpsq4,'k');
    t = linspace(3e-2,4.2e-1,2);
    loglog(t,t*1.2,'--k','LineWidth',.4);
    hold off;

    axis([2e-2 2 1e-6 1e0]);
    ylabel('$\langle |\delta\!B_{\!\perp}|^{2} \rangle$');
    set(gca,'XTick',[1e-1 1e0],'XTickLabel',[ ]);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'Units','normalized','Position',[0.15,0.54,0.84,0.43]);
    plotTickLatex2D('ytickdx',-0.007,'ylabeldx',-0.025,'ytickdy',0.004);

    ah1 = gca;
    leg1 = legend(ah1,h(1:4),'$~S = 0.003$','$~S = 0.001$', ...
        '$~S = 0.0003$','$~S = 0.0001$','Location',[0.635 0.57 0.32 0.1783]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');

    text(7e-2,2.6e-1,'$\propto\! t$','FontSize',10,'Interpreter','latex');
    
    
    subplot(2,1,2);
    
    loglog(time1,dBprlsq1,'Color',[232/255 0 0]);
    hold on;
    loglog(time2,dBprlsq2,'Color',[34/255 139/255 34/255]);
    loglog(time3,dBprlsq3,'b');
    loglog(time4,dBprlsq4,'k');
    %loglog(time5,dBprlsq5,'Color',[1 128/255 40/255]);
    
    t = linspace(3e-2,5e-1,2);
    loglog(t,t.^2*0.3,'--k','LineWidth',.4);
    hold off;

    axis([2e-2 2 1e-6 1e0]);
    xlabel('$S t$'); ylabel('$\langle \delta\!B^2_{||}\, \rangle$');
    set(gca,'XTick',[1e-1 1e0]);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'Units','normalized','Position',[0.15,0.09,0.84,0.43]);
    plotTickLatex2D('xlabeldy',0.015,'ytickdx',-0.007,'ylabeldx',-0.018,'ytickdy',0.004);

    text(7e-2,7e-3,'$\propto\! t^2$','FontSize',10,'Interpreter','latex');
    
    %}
    
   
    
    fig(33); clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,5.5]);
    set(gcf,'PaperPositionMode','auto');
    
    S = [1,3,10,30]*1e-4; SD = [0.02886,0.06423,0.1412,0.2784];
  
    plot(time1,-Lam1,'Color',[232/255 0 0]);
    hold on;
    plot(0.249,SD(4),'+','Color',[232/255 0 0]);
    plot(time2,-Lam2,'Color',[34/255 139/255 34/255]);
    plot(0.122,SD(3),'+','Color',[34/255 139/255 34/255]);
    plot(time3,-Lam3,'b');
    plot(0.0573,SD(2),'+b');
    plot(time4,-Lam4,'k');
    plot(0.0296,SD(1),'+k');
    %plot(time5,Lam5,'Color',[1 128/255 40/255]);
    hold off;
    axis([0 1 -0.06 0.3]);
    %axis([0 1 -0.068 0.008]);
    set(gca,'TickLength',[0.025 0.025],'Fontsize',10);
    set(gca,'YTick',[-0.1 0 0.1 0.2 0.3]); 
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.82]);
    xlabel('$S t$'); ylabel('$\langle\, \Lambda_{\rm f} \, \rangle$');%\langle\, \beta_\perp / \beta_{||}\, + 2/\beta_{||}\, - 1 \rangle$')
    plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.004,'ylabeldx',-0.045,'ytickdy',0.01);
    set(gca,'XMinorTick','on','YMinorTick','on')
    text(0.043,0.26,'(b)','Fontsize',10);
    
    %{
    hax = axes('Position',[0.59 0.47 0.33 0.4]);
    SL = linspace(-4.2,-2.3,2);
    plot(SL,2/3*SL+1.14,'--k','LineWidth',0.5);
    hold on;
    plot(log10(S(1)),log10(SD(1)),'+k');
    plot(log10(S(2)),log10(SD(2)),'+b');
    plot(log10(S(3)),log10(SD(3)),'+','Color',[34/255 139/255 34/255]);
    plot(log10(S(4)),log10(SD(4)),'+','Color',[232/255 0 0]);
    
    axis([-4.5 -2 -1.8 -0.2]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10);
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel('${\rm lg}~S$'); ylabel('${\rm lg}~\langle\,\Lambda_{\rm f}\,\rangle_{\rm max}$');
    plotTickLatex2D('xlabeldy',-0.015,'xtickdy',-0.01,'ytickdx',-0.005,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(-3.85,-0.7,'$\propto\!S^{2/3}$','FontSize',10,'Interpreter','latex');    
    hold off;
    %}
    
    hax = axes('Position',[0.59 0.47 0.33 0.4]);
    SL = linspace(-4.2,-2.3,2);
    loglog(10.^SL,10.^(2/3*SL+1.14),'--k','LineWidth',0.5);
    hold on;
    loglog((S(1)),(SD(1)),'+k');
    loglog((S(2)),(SD(2)),'+b');
    loglog((S(3)),(SD(3)),'+','Color',[34/255 139/255 34/255]);
    loglog((S(4)),(SD(4)),'+','Color',[232/255 0 0]);
    
    axis([3e-5 1e-2 0.01 0.8]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10,'XTick',[1e-4 1e-3 1e-2]);
    set(gca,'XMinorTick','on','YMinorTick','on','YTick',[1e-2 1e-1]);
    xlabel('$S$'); ylabel('$\langle\,\Lambda_{\rm f}\,\rangle_{\rm max}$');
    plotTickLatex2D('xlabeldy',-0.015,'xtickdy',-0.01,'ytickdx',-0.002,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(10^-3.85,10^-0.65,'$\propto\!S^{2/3}$','FontSize',10,'Interpreter','latex');    
    hold off;
    
    
    % scattering plots
    fig(333);clf;
    
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,5.5]);
    set(gcf,'PaperPositionMode','auto');
    
    for i=1:length(time1)-1
        dlnBdt1(i) = 0.5*(log(ME1(i+1))-log(ME1(i)));
        dlnDdt1(i) = (Del1(i+1)/Del1(i)-1);
    end
    dlnBdt1(i+1) = dlnBdt1(i);
    dlnDdt1(i+1) = dlnDdt1(i);
    for i=1:length(time2)-1
        dlnBdt2(i) = 0.5*(log(ME2(i+1))-log(ME2(i)));
        dlnDdt2(i) = (Del2(i+1)/Del2(i)-1);
    end
    dlnBdt2(i+1) = dlnBdt2(i);
    dlnDdt2(i+1) = dlnDdt2(i);
    for i=1:length(time3)-1
        dlnBdt3(i) = 0.5*(log(ME3(i+1))-log(ME3(i)));
        dlnDdt3(i) = (Del3(i+1)/Del3(i)-1);
    end
    dlnBdt3(i+1) = dlnBdt3(i);
    dlnDdt3(i+1) = dlnDdt3(i);
    for i=1:length(time4)-1
        dlnBdt4(i) = 0.5*(log(ME4(i+1))-log(ME4(i)));
        dlnDdt4(i) = (Del4(i+1)/Del4(i)-1);
    end
    dlnBdt4(i+1) = dlnBdt4(i);
    dlnDdt4(i+1) = dlnDdt4(i);
    %{
    for i=1:length(time5)-1
        dlnBdt5(i) = 0.5*(log(ME5(i+1))-log(ME5(i)));
        dlnDdt5(i) = (Del5(i+1)/Del5(i)-1);
    end
    dlnBdt5(i+1) = dlnBdt5(i);
    dlnDdt5(i+1) = dlnDdt5(i);
    %}
    
    nu1 = ( 3*dlnBdt1'./Del1 - dlnDdt1' );
    nu2 = ( 3*dlnBdt2'./Del2 - dlnDdt2' );
    nu3 = ( 3*dlnBdt3'./Del3 - dlnDdt3' );
    nu4 = ( 3*dlnBdt4'./Del4 - dlnDdt4' );
    %nu5 = ( 3*dlnBdt5'./Del5 - dlnDdt5' );
    
    %{
    notdone=1; i=length(dlnBdt1);
    while notdone
        if (dlnBdt1(i) > 0)
            notdone = 0;
        else
            i=i-1;
        end
    end
    nu1(1:i) = -1;
    notdone=1; i=length(dlnBdt2);
    while notdone
        if (dlnBdt2(i) > 0)
            notdone = 0;
        else
            i=i-1;
        end
    end
    nu2(1:i) = -1;
    notdone=1; i=length(dlnBdt3);
    while notdone
        if (dlnBdt3(i) > 0)
            notdone = 0;
        else
            i=i-1;
        end
    end
    nu3(1:i) = -1;
    notdone=1; i=length(dlnBdt4);
    while notdone
        if (dlnBdt4(i) > 0)
            notdone = 0;
        else
            i=i-1;
        end
    end
    nu4(1:i) = -1; 
    %}
    
    plot(time1,nu1,'Color',[232/255 0 0]);
    hold on;
    plot(time2,nu2,'Color',[34/255 139/255 34/255]);
    plot(time3,nu3,'b');
    plot(time4,nu4,'k');
    %plot(time5,nu5,'c');
    hold off;
    axis([0 1.2 0 0.4]);
    set(gca,'TickLength',[0.025 0.025],'Fontsize',10);
    set(gca,'XTick',[0.2 0.4 0.6 0.8 1.0],'XMinorTick','on','YMinorTick','on');
    set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.82]);
    xlabel('$S t$'); ylabel('$\nu_{\rm eff}~/ S \beta_0$');
    plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.004,'ylabeldx',-0.045,'ytickdy',0.01);
    
    %{
    plot(time1,nu1./num1,'Color',[232/255 0 0]);
    hold on;
    plot(time2,nu2./num2,'Color',[34/255 139/255 34/255]);
    plot(time3,nu3./num3,'b');
    plot(time4,nu4./num4,'k');
    hold off;
    axis([0 2 1e-3 3]);
    set(gca,'TickLength',[0.025 0.025],'Fontsize',10);
    %set(gca,'YTick',[0 0.1 0.2],'XMinorTick','on','YMinorTick','on');
    set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.82]);
    xlabel('$S t$'); ylabel('$\nu_{\rm eff}$');
    plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.004,'ylabeldx',-0.045,'ytickdy',0.01);
%}
    
%{
    ah1 = gca;
    leg1 = legend(ah1,'$~S = 0.003$','$~S = 0.001$', ...
        '$~S = 0.0003$','$~S = 0.0001$','Location','NorthWest');%[0.635 0.57 0.32 0.1783]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
%}
        
        
elseif i==4
%
% MIRROR
% box average evolution of magnetic energy; Reynolds/Maxwell/pressure 
% stresses; instability parameter
%
 

    shear = 0.003;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-fast/';        % directory
    fname = 'mrshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam1 = C{21};
    Del1 = C{22};
    
    time1 = time*shear;

    ME1 = MEx+MEy+MEz;
    B1sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time1).^2);
    by = (pitch-time1)./sqrt(1.0+(pitch-time1).^2);
    dBprlsq1 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq1 = 2*ME1-B1sq-dBprlsq1;
        
    beta = 100./ME1;   
    for i=1:length(time1)-1
        dlnBdt1(i) = 0.5*(log(ME1(i+1))-log(ME1(i)));
    end
    dlnBdt1(i+1) = dlnBdt1(i);
    num1 = 3*beta(:).*dlnBdt1(:); %shear*beta.*bx.*by;
    
    shear = 1e-3;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-medium/';        % directory
    fname = 'mrshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam2 = C{21};
    Del2 = C{22};
    
    time2 = time*shear;

    ME2 = MEx+MEy+MEz;
    B2sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time2).^2);
    by = (pitch-time2)./sqrt(1.0+(pitch-time2).^2);
    dBprlsq2 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq2 = 2*ME2-B2sq-dBprlsq2;
   
    
    beta = 100./ME2;   
    for i=1:length(time2)-1
        dlnBdt2(i) = 0.5*(log(ME2(i+1))-log(ME2(i)));
    end
    dlnBdt2(i+1) = dlnBdt2(i);
    num2 = 3*beta(:).*dlnBdt2(:); %shear*beta.*bx.*by;
    
    
    shear = 3e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/';        % directory
    fname = 'mrshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam3 = C{21};
    Del3 = C{22};
    
    time3 = time*shear;

    ME3 = MEx+MEy+MEz;
    B3sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time3).^2);
    by = (pitch-time3)./sqrt(1.0+(pitch-time3).^2);
    dBprlsq3 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq3 = 2*ME3-B3sq-dBprlsq3;
    
    
    beta = 100./ME3;   
    for i=1:length(time3)-1
        dlnBdt3(i) = 0.5*(log(ME3(i+1))-log(ME3(i)));
    end
    dlnBdt3(i+1) = dlnBdt3(i);
    num3 = 3*beta(:).*dlnBdt3(:); %shear*beta.*bx.*by;
    
    
    shear = 1e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-super-slow/';        % directory
    fname = 'mrshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    MEx  = C{11};
    MEy  = C{12};
    MEz  = C{13};
    Bxa  = C{17};
    Bya  = C{18};
    Maxw = C{20};
    Lam4 = C{21};
    Del4 = C{22};
    
    time4 = time*shear;

    ME4 = MEx+MEy+MEz;
    B4sq = Bxa.^2+Bya.^2;
    
    pitch = Bya(1)/Bxa(1);
    bx = 1.0./sqrt(1.0+(pitch-time4).^2);
    by = (pitch-time4)./sqrt(1.0+(pitch-time4).^2);
    dBprlsq4 = bx.^2.*(2*MEx-Bxa.^2)+by.^2.*(2*MEy-Bya.^2)-2*bx.*by.*(Maxw+Bxa.*Bya);
    dBprpsq4 = 2*ME4-B4sq-dBprlsq4;
      
    beta = 100./ME4;   
    for i=1:length(time4)-1
        dlnBdt4(i) = 0.5*(log(ME4(i+1))-log(ME4(i)));
    end
    dlnBdt4(i+1) = dlnBdt4(i);
    num4 = 3*beta(:).*dlnBdt4(:); %shear*beta.*bx.*by;
    
    
    %
    % if you want one panel with dBprl
    %
    figure(4);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    
    set(gcf,'Position',[1,1,11.5,5]);
%    set(gcf,'Position',[1,1,11.5,5.5]);

    S = [1e-4,3e-4,1e-3,3e-3]; SD = [0.1382,0.2175,0.35,0.555];
    BL = linspace(-4.2,-2.55,2); BD = [0.04201,0.07083,0.1301,0.1439];

    h(1) = loglog(time1,dBprlsq1,'Color',[232/255 0 0]);
    hold on;
    %loglog(SD(4),BD(4),'+','Color',[232/255 0 0]);
    h(2) = loglog(time2,dBprlsq2,'Color',[34/255 139/255 34/255]);
    %loglog(SD(3),BD(3),'+','Color',[34/255 139/255 34/255]);
    h(3) = loglog(time3,dBprlsq3,'b');
    %loglog(SD(2),BD(2),'+b');
    h(4) = loglog(time4,dBprlsq4,'k');
    %loglog(SD(1),BD(1),'+k');
    t = linspace(1e-1,1.05,2);
    loglog(t,t.^(4/3)*1.2,'--k','LineWidth',.5);
    hold off;
    text(7.1e-2,4e-1,'(a)','Fontsize',10);
    axis([6e-2 3 1e-5 2]);
    xlabel('$S t$'); ylabel('$\langle \delta\!B^2_{||}\, \rangle / B^2_0$');
    set(gca,'XTick',[1e-1 1e0]);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'Units','normalized','Position',[0.15 0.145 0.84 0.83]);
    %set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.83]);
    plotTickLatex2D('xtickdy',-0.015,'xlabeldy',-0.02,'ytickdx',-0.007,'ylabeldx',-0.018,'ytickdy',0.01);
    %plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.007,'ylabeldx',-0.018,'ytickdy',0.01);
    text(2e-1,0.6,'$\propto\! t^{4/3}$','FontSize',10,'Interpreter','latex');
    
    ah1 = gca;
    leg1 = legend(ah1,h(1:4),'$~S = 0.003$','$~S = 0.001$', ...
        '$~S = 0.0003$','$~S = 0.0001$','Location',[0.635 0.3 0.32 0.1783]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');

    
    %{
    hax = axes('Position',[0.65 0.34 0.28 0.37]);
    plot(BL,0.5*BL+0.615,'--k','LineWidth',0.5);
    hold on;
    plot(log10(S(1)),log10(BD(1)),'+k');
    plot(log10(S(2)),log10(BD(2)),'+b');
    plot(log10(S(3)),log10(BD(3)),'+','Color',[34/255 139/255 34/255]);
    plot(log10(S(4)),log10(BD(4)),'+','Color',[232/255 0 0]);
    
    axis([-4.5 -2 -1.6 -0.6]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10);
    set(gca,'XMinorTick','on','YMinorTick','on','YTick',[-1.5 -1]);
    xlabel('${\rm lg}~S$'); title('${\rm lg}~\langle \delta\!B^2_{||}\, \rangle_{\rm lin,sat}$');
    plotTickLatex2D('xlabeldy',-0.01,'xtickdy',-0.01,'ytickdx',-0.005,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(-4.1,-0.85,'$\propto\!S^{1/2}$','FontSize',10,'Interpreter','latex');    
    hold off;
    %}
    
    
    fig(44); clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Position',[1,1,11.5,5]);
    %set(gcf,'Position',[1,1,11.5,5.5]);
    
    S = [1,3,10,30]*1e-4; SD = [0.1041,0.1791,0.3254,0.5677];
    
    plot(time1,Lam1,'Color',[232/255 0 0]);
    hold on;
    plot(0.352,SD(4),'+','Color',[232/255 0 0]);
    plot(time2,Lam2,'Color',[34/255 139/255 34/255]);
    plot(0.233,SD(3),'+','Color',[34/255 139/255 34/255]);
    plot(time3,Lam3,'b');
    plot(0.144,SD(2),'+b');
    plot(time4,Lam4,'k');
    plot(0.0933,SD(1),'+k');
    hold off;
    text(0.12,0.51,'(b)','Fontsize',10);
    axis([0 3 -0.05 0.6]);
    set(gca,'TickLength',[0.025 0.025],'Fontsize',10);
    set(gca,'YTick',[ 0 0.1 0.2 0.3 0.4 0.5 0.6]); set(gca,'XTick',[0 0.5 1 1.5 2 2.5 3]);
    set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.83]);
    %set(gca,'Units','normalized','Position',[0.15,0.14,0.84,0.82]);
    xlabel('$S t$'); ylabel('$\langle\, \Lambda_{\rm m} \, \rangle$');
    plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.004,'ylabeldx',-0.046,'ytickdy',0.01);
    set(gca,'XMinorTick','on','YMinorTick','on');
    
    %{
    hax = axes('Position',[0.59 0.427 0.33 0.45]);
    %hax = axes('Position',[0.59 0.47 0.33 0.4]);
    SL = linspace(-4.14,-2.3,2);
    plot(SL,0.5*SL+1.015,'--k','LineWidth',0.5);
    hold on;
    plot(log10(S(1)),log10(SD(1)),'+k');
    plot(log10(S(2)),log10(SD(2)),'+b');
    plot(log10(S(3)),log10(SD(3)),'+','Color',[34/255 139/255 34/255]);
    plot(log10(S(4)),log10(SD(4)),'+','Color',[232/255 0 0]);
    
    axis([-4.5 -2 -1.2 0]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10);
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel('${\rm lg}~S$'); ylabel('${\rm lg}~\langle\,\Lambda_{\rm m}\,\rangle_{\rm max}$');
    plotTickLatex2D('xlabeldy',-0.015,'xtickdy',-0.01,'ytickdx',-0.005,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(-3.85,-0.4,'$\propto\!S^{1/2}$','FontSize',10,'Interpreter','latex');    
    hold off;
    %}
    
    hax = axes('Position',[0.59 0.427 0.33 0.45]);
    %hax = axes('Position',[0.59 0.47 0.33 0.4]);
    SL = linspace(-4.14,-2.3,2);
    loglog(10.^SL,10.^(0.5*SL+1.015),'--k','LineWidth',0.5);
    hold on;
    loglog((S(1)),(SD(1)),'+k');
    loglog((S(2)),(SD(2)),'+b');
    loglog((S(3)),(SD(3)),'+','Color',[34/255 139/255 34/255]);
    loglog((S(4)),(SD(4)),'+','Color',[232/255 0 0]);
    
    axis([3e-5 1e-2 0.06 1]);
    set(gca,'TickLength',[0.04 0.04],'FontSize',10);
    set(gca,'XMinorTick','on','YMinorTick','on','XTick',[1e-4 1e-3 1e-2]);
    xlabel('$S$'); ylabel('$\langle\,\Lambda_{\rm m}\,\rangle_{\rm max}$');
    plotTickLatex2D('xlabeldy',-0.015,'xtickdy',-0.01,'ytickdx',-0.002,'ylabeldx',-0.004,'xtickdx',-0.005,'ytickdy',0.005);
    text(10^-3.85,10^-0.4,'$\propto\!S^{1/2}$','FontSize',10,'Interpreter','latex');    
    hold off;



%{
    figure(14);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,5.5]);
    set(gcf,'PaperPositionMode','auto');
    
    shear = 3e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/';        % directory
    fname = 'vtk/mrshear';
    
    fn1=0;
    fn2=628;
    Nx = 1152; Ny = 1152;

    for ff=fn1:fn2
        
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

        [time,var,name,status] = readvtk(Grid,filename,2);
        U.d = squeeze(var);

        cgl = (U.bx.^2 + U.by.^2 + U.bz.^2) ./ (U.d.^(4/3));
        cgl = sum(sum(cgl))/Nx/Ny;
        cgl = log(cgl);
        
        c(ff+1) = cgl;
        
    end
    
    for i=fn1+1:fn2 
        change(i+1) = c(i+1)-c(i);
    end
    
    plot(change);
    %}

    
elseif i==5
%
%
%    

Lx = 1152; Ly = 1152; shear = 3e-4;
dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/vtk/';        % directory
fname = 'mrshear';
    
f = 50;

    figure(5);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,12.6]);
    set(gcf,'PaperPositionMode','auto');
           
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
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
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = -1/sqrt(5);
    bx  = 2/sqrt(5);
    by  = by0 - bx*shear*time;
    dBy = U.by-by;
    dBx = U.bx-bx;
    dBprl = (by*dBy+bx*dBx)/sqrt(by*by+bx*bx);
    
    %{
    dBprlp = remap2d(Grid.x1,Grid.x2,shear,time,-1,dBprl);
    dBk = 2*abs(fft2(dBprlp)/Grid.nx2/Grid.nx1).^2;
    kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
    kx = kx1 * (2*pi/Grid.dx1);
    ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1 ) - 1/2;
    ky = ky1 * (2*pi/Grid.dx2);
    [KX,KY] = meshgrid(kx,ky);
    tremap = mod(time + Ly / (2.0 * shear * Lx) , Ly / (shear * Lx)) - Ly / (2.0 * shear * Lx);
    KX = KX + KY*shear*tremap;
    Kprl = abs(bx*KX+by*KY)/sqrt(bx*bx+by*by);
    Kprp = abs(by*KX-bx*KY)/sqrt(bx*bx+by*by);
    nbins = 576;
    bins = max(max(max(Kprl)),max(max(Kprp)))*linspace(0,1,nbins);
    bins = log10(bins+2*pi/Lx/sqrt(2));
    [bincounts,ind]= histc(log10(Kprl),bins);
    dBkprl = zeros(nbins,1);
    for i=1:Grid.nx1
      for j=1:Grid.nx2
        if ind(i,j) ~= 0
            dBkprl(ind(i,j)) = dBkprl(ind(i,j)) + dBk(i,j);
        end
      end
    end
    dBkprl = 0.5*dBkprl;
    [bincounts,ind]= histc(log10(Kprp),bins);
    dBkprp = zeros(nbins,1);
    for i=1:Grid.nx1
      for j=1:Grid.nx2
          if ind(i,j) ~= 0
            dBkprp(ind(i,j)) = dBkprp(ind(i,j)) + dBk(i,j);
          end
      end
    end
    dBkprp = 0.5*dBkprp;
    kprl = 10.^bins;%'*sqrt(200);
    loglog(kprl,kprl.^2.*dBkprl','k');
    hold on;
    kprp = 10.^bins;%'*sqrt(200);
    loglog(kprp,kprp.^2.*dBkprp','r');
    loglog(kprp,3e-7*kprp.^(-8/3),'--b');
    %loglog(kprl,3e-1*kprl.^(-4),'b');
    loglog(1/sqrt(200)*[1 1],[1e-8 1],'--k','LineWidth',0.5);
    loglog([1 1],[1e-8 1],'--k','LineWidth',0.5);
    hold off;
    %}
    
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx1 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx1 = kx(1:nsx/2) * (2*pi/pstep);

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx1;			% step taken in y direction during field-line trace
    nsy =  2*Ly/Grid.dx1;			% number of steps taken to get periodic series
    nsx =   Lx/Grid.dx1;			% number of steps taken to shift x0 and cover space
    zz = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end
    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy1 = sum(zfy,2)/nsx/nsx;
    ky1 = ky*2*pi;
   
    
    
    f = 333;
    
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
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
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = -1/sqrt(5);
    bx  = 2/sqrt(5);
    by  = by0 - bx*shear*time;
    dBy = U.by-by;
    dBx = U.bx-bx;
    dBprl = (by*dBy+bx*dBx)/sqrt(by*by+bx*bx);
    
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx2 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx2 = kx(1:nsx/2) * (2*pi/pstep);

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx1;			% step taken in y direction during field-line trace
    nsy =  2*Ly/Grid.dx1;			% number of steps taken to get periodic series
    nsx =  Lx/Grid.dx1;			% number of steps taken to shift x0 and cover space
    zz = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    end

    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy2 = sum(zfy,2)/nsx/nsx;
    ky2 = ky*2*pi;
    
    
    f = 628;
    
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
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
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = -1/sqrt(5);
    bx  = 2/sqrt(5);
    by  = by0 - bx*shear*time;
    dBy = U.by-by;
    dBx = U.bx-bx;
    dBprl = (by*dBy+bx*dBx)/sqrt(by*by+bx*bx);
    
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx3 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx3 = kx(1:nsx/2) * (2*pi/pstep);

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx1;			% step taken in y direction during field-line trace
    nsy =  2*Ly/Grid.dx1;			% number of steps taken to get periodic series
    nsx =  Lx/Grid.dx1;			% number of steps taken to shift x0 and cover space
    zz = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    end

    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy3 = sum(zfy,2)/nsx/nsx;
    ky3 = ky*2*pi;
    
    
    subplot(3,1,1);
    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.073*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    axis([2e-3 2e0 1e-6 5e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XTickLabel',['']);
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.735,0.84,0.26]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.0);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$~k_{||}$','$~k_\perp$','$~\rho^{-1}_{\rm i}$','Location','Northeast');%[0.77 0.88 0.1928 0.096]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(3.1e-3,5e-2,'$St = 0.15$','FontSize',10,'Interpreter','latex');
    text(0.028,4.2e-2,'$k^{\rm max}_{||}$','FontSize',10,'Interpreter','latex');
    text(0.0573,4.9e-2,'$k^{\rm max}_{\perp}$','FontSize',10,'Interpreter','latex');
    arrow([0.036 1.2e-2],[0.036 1.8e-3],5,'BaseAngle',60,'EdgeColor','b','FaceColor','b','LineWidth',1.0);
    arrow([0.064 1.2e-2],[0.064 1.8e-3],5,'BaseAngle',60,'EdgeColor',[232/255 0 0],'FaceColor',[232/255 0 0],'LineWidth',1.0);

    subplot(3,1,2);
    
    h(4) = loglog(kx2,zfx2,'b');
    hold on;
    h(5) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    loglog(0.107*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    kk = linspace(3.5e-2,1.07e-1,2);
    loglog(kk,3.5e-8*kk.^(-11/3),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,2.8e-9*kk.^(-4.8),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,4e-1,2);
    loglog(kk,2.2e-7*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-7 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XTickLabel',[ ]);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.401,0.84,0.319]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.0);
    text(1.8e-2,3e-4,'$\propto\! k^{-11/3}$','FontSize',10,'Interpreter','latex');
    text(4e-2,9e-6,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(2.2e-1,1.6e-3,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(3.1e-3,8e-7,'$St = 1$','FontSize',10,'Interpreter','latex');
    
    subplot(3,1,3);
    
    h(6) = loglog(kx3,zfx3,'b');
    hold on;
    loglog(0.16*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    h(7) = loglog(ky3,zfy3,'Color',[232/255 0 0]); 
    kk = linspace(4.2e-2,1.6e-1,2);
    loglog(kk,5.5e-8*kk.^(-11/3),'k','LineWidth',0.6);
    kk = linspace(1.7e-1,2.8e-1,2);
    loglog(kk,4.2e-9*kk.^(-4.8),'k','LineWidth',0.6);
    kk = linspace(1.7e-1,4.4e-1,2);
    loglog(kk,6e-7*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-7 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    xlabel('$k$');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.067,0.84,0.319]);
    plotTickLatex2D('xtickdy',0.008,'xlabeldy',0.028,'ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.0);
    text(6e-2,4e-6,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(2.7e-1,1.6e-3,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(2e-2,2e-4,'$\propto\! k^{-11/3}$','FontSize',10,'Interpreter','latex');
    text(3.1e-3,8e-7,'$St = 1.9$','FontSize',10,'Interpreter','latex');
        

    figure(55);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,8.35]);
    set(gcf,'PaperPositionMode','auto');
    
    subplot(2,1,1);
    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.073*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    axis([2e-3 2e0 1e-6 5e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XTickLabel',['']);
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.602,0.84,0.393]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$~k_{||}$','$~k_\perp$','$~\rho^{-1}_{\rm i}$','Location','Northeast');%[0.77 0.82 0.194 0.116]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(3.1e-3,5e-2,'$St = 0.15$','FontSize',10,'Interpreter','latex');
    text(0.028,4.2e-2,'$k^{\rm max}_{||}$','FontSize',10,'Interpreter','latex');
    text(0.0573,4.9e-2,'$k^{\rm max}_{\perp}$','FontSize',10,'Interpreter','latex');
    arrow([0.036 1.2e-2],[0.036 1.8e-3],5,'BaseAngle',60,'EdgeColor','b','FaceColor','b','LineWidth',1.0);
    arrow([0.064 1.2e-2],[0.064 1.8e-3],5,'BaseAngle',60,'EdgeColor',[232/255 0 0],'FaceColor',[232/255 0 0],'LineWidth',1.0);
    
    subplot(2,1,2);
    
    h(4) = loglog(kx2,zfx2,'b');
    hold on;
    h(5) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    loglog(0.107*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    kk = linspace(3.5e-2,1.07e-1,2);
    loglog(kk,3.2e-8*kk.^(-11/3),'k','LineWidth',0.4);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,2.5e-9*kk.^(-4.8),'k','LineWidth',0.4);
    kk = linspace(1.1e-1,4e-1,2);
    loglog(kk,2.2e-7*kk.^(-4.8),'k','LineWidth',0.4);
    hold off;
    axis([2e-3 2e0 1e-7 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    xlabel('$k$');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.095,0.84,0.481]);
    plotTickLatex2D('xtickdy',-0.0,'xlabeldy',0.01,'ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    text(1.8e-2,3e-4,'$\propto\! k^{-11/3}$','FontSize',10,'Interpreter','latex');
    text(4e-2,9e-6,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(2.2e-1,1.6e-3,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(3.1e-3,8e-7,'$St = 1$','FontSize',10,'Interpreter','latex');
    
        
    
    figure(555);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,7.7]);
    set(gcf,'PaperPositionMode','auto');
    
    subplot(2,1,1);
    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.073*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    axis([2e-3 2e0 1e-6 5e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XTickLabel',['']);
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.572,0.84,0.425]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$~k_{||}$','$~k_\perp$','$~\rho^{-1}_{\rm i}$','Location','Northeast');%[0.77 0.82 0.194 0.116]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(3.1e-3,5e-2,'$St = 0.15$','FontSize',10,'Interpreter','latex');
    text(0.028,4.2e-2,'$k^{\rm max}_{||}$','FontSize',10,'Interpreter','latex');
    text(0.0573,4.9e-2,'$k^{\rm max}_{\perp}$','FontSize',10,'Interpreter','latex');
    arrow([0.036 1.2e-2],[0.036 1.8e-3],5,'BaseAngle',60,'EdgeColor','b','FaceColor','b','LineWidth',1.0);
    arrow([0.064 1.2e-2],[0.064 1.8e-3],5,'BaseAngle',60,'EdgeColor',[232/255 0 0],'FaceColor',[232/255 0 0],'LineWidth',1.0);
    
    subplot(2,1,2);
    
    h(4) = loglog(kx2,zfx2,'b');
    hold on;
    h(5) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    loglog(0.107*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    kk = linspace(3.5e-2,1.07e-1,2);
    loglog(kk,3.2e-8*kk.^(-11/3),'k','LineWidth',0.4);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,2.5e-9*kk.^(-4.8),'k','LineWidth',0.4);
    kk = linspace(1.1e-1,4e-1,2);
    loglog(kk,2.2e-7*kk.^(-4.8),'k','LineWidth',0.4);
    hold off;
    axis([2e-3 2e0 1e-6 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    xlabel('$k$');
    ylabel('$|\delta\!B_{||,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.099,0.84,0.452]);
    plotTickLatex2D('xtickdy',-0.0,'xlabeldy',0.005,'ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    text(1.8e-2,3e-4,'$\propto\! k^{-11/3}$','FontSize',10,'Interpreter','latex');
    text(4e-2,9e-6,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(2.2e-1,1.6e-3,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(3.1e-3,8e-6,'$St = 1$','FontSize',10,'Interpreter','latex');
    
        
    
   
elseif i==6
    
Lx = 1152; Ly = 1152; shear = 3e-4;

%dir ='/Volumes/My Passport Studio/pegasus/fhs-super-slow/vtk/';
dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/vtk/';        % directory
fname = 'fhshear';

    f=22;
    
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
    end

    % declare file name
    filename = [dir,fname,'.',numlab,'.vtk'];

    % open file and initialize grid
    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);

    [time,var,name,status] = readvtk(Grid,filename,1);
    U.bz  = squeeze(var(3,:,:,:));
    
    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d = squeeze(var)-1;
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = 3/sqrt(13);
    bx  = 2/sqrt(13);
    by  = by0 - bx*shear*time;
    dBz = U.bz;
    
    %{
    dBzp = remap2d(Grid.x1,Grid.x2,shear,time,-1,dBz);
    dBk = 2*abs(fft2(dBzp)/Grid.nx1/Grid.nx2).^2;
    kx1 = mod( 1/2 + (0:(Grid.nx1-1))/Grid.nx1 , 1 ) - 1/2;
    kx = kx1 * (2*pi/Grid.dx1);
    ky1 = mod( 1/2 + (0:(Grid.nx2-1))/Grid.nx2 , 1) - 1/2;
    ky = ky1 * (2*pi/Grid.dx2);
    [KX,KY] = meshgrid(kx,ky);
    KX = KX + KY*shear*tremap;
    Kprl = abs(bx*KX+by*KY)/sqrt(bx*bx+by*by);
    Kprp = abs(by*KX-bx*KY)/sqrt(bx*bx+by*by);
    KK   = sqrt(KX.^2+KY.^2);
    nbins = 576;
    bins = max(max(max(Kprl)),max(max(Kprp)))*linspace(0,1,nbins);
    bins = log10(bins+2*pi/Lx/sqrt(2));
    [bincounts,ind]= histc(log10(Kprl),bins);
    dBkprl = zeros(nbins,1);
    for i=1:Grid.nx1
      for j=1:Grid.nx2
        if ind(i,j) ~= 0
            dBkprl(ind(i,j)) = dBkprl(ind(i,j)) + dBk(i,j);
        end
      end
    end
    dBkprl = dBkprl;
    [bincounts,ind]= histc(log10(Kprp),bins);
    dBkprp = zeros(nbins,1);
    for i=1:Grid.nx1
      for j=1:Grid.nx2
          if ind(i,j) ~= 0
            dBkprp(ind(i,j)) = dBkprp(ind(i,j)) + dBk(i,j);
          end
      end
    end
    kprl = 10.^bins;%'*sqrt(200);
    loglog(kprl,dBkprl','k');
    hold on;
    kprp = 10.^bins;%'*sqrt(200);
    loglog(kprp,dBkprp','r');
    %loglog(kprl,dBkk,'b');
    %loglog(kprp,6e-8*kprp.^(-2.8),'--b');
    %loglog(kprp,2e8*kprp.^(-3),'b');
    loglog(1/sqrt(200)*[1 1],[1e-8 1],'--k','LineWidth',0.5);
    loglog([1 1],[1e-8 1],'--k','LineWidth',0.5);
    hold off;
    %}
            
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx1 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx1 = kx(1:nsx/2) * (2*pi/pstep);

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx2;			% step taken in y direction during field-line trace
    nsy = 2*Ly/Grid.dx1;			% number of steps taken to get periodic series
    nsx =   Lx/Grid.dx2;			% number of steps taken to shift x0 and cover space
    zz = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    %plot(xi,yi,'.',[x0],[y0],'.r');
    %axis([-580 580 -580 580]);
    %pause(0.5);

    end
    
    %{
    zfy = zeros(nsy,1);
    for n=1:nsx
        zzf = fft(zz(:,n))/nsy/2;
        zfy = zfy + 2*abs(zzf(1:nsy)).^2;
    end
    zfy1 = zfy/nsx;

    ky = mod( 1/2 + (0:(2*nsy-1))/2/nsy , 1 ) - 1/2;
    ky1 = ky(1:nsy) * (2*pi/pstep);
    %}
    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy1 = sum(zfy,2)/nsx/nsx;
    ky1 = ky*2*pi;
    
    
    f=333;
        
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
    end

    % declare file name
    filename = [dir,fname,'.',numlab,'.vtk'];

    % open file and initialize grid
    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);

    [time,var,name,status] = readvtk(Grid,filename,1);
    U.bz  = squeeze(var(3,:,:,:));
    
    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d   = squeeze(var);
    dn = U.d - 1;
    
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = 3/sqrt(13);
    bx  = 2/sqrt(13);
    by  = by0 - bx*shear*time;
    dBz = U.bz;
           
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx2 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx2 = kx(1:nsx/2) * (2*pi/pstep);

% Calculate the frequencies

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx2;			% step taken in y direction during field-line trace
    nsy = 2*Ly/Grid.dx2;			% number of steps taken to get periodic series
    nsx = Lx/Grid.dx1;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsy,nsx);
    

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));
    
    ni = interp2(Grid.x1,Grid.x2,dn',xi,yi);
    nn(1:nsy,n) = inpaint_nans(ni(1:nsy));

    end

    %{
    zfy = zeros(nsy/2,1);
    for n=1:nsx
        zzf = fft(zz(:,n))/nsy;
        zfy = zfy + 2*abs(zzf(1:nsy/2)).^2;
    end
    zfy2 = zfy/nsx;

    ky = mod( 1/2 + (0:(nsy-1))/nsy , 1 ) - 1/2;
    ky2 = ky(1:nsy/2) * (2*pi/pstep);
    %}
    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy2 = sum(zfy,2)/nsx/nsx;
    ky2 = ky*2*pi;
    
    
    figure(6);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,8.2]);
    set(gcf,'PaperPositionMode','auto');
    
    subplot(2,1,1);

    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.0696*[1 1],[1e-8 1e-1],'--k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-7 1e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XTickLabel',['']);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'YMinorTick','on');
    ylabel('$|\delta\!B_{z\!,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.54,0.84,0.43]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$~k_{||}$','$~k_\perp$','$~\rho_{\rm i}^{-1}$','Location','Northeast');%[0.77 0.807 0.1934 0.116]);
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(3.1e-3,1e-2,'$St = 0.066$','FontSize',10,'Interpreter','latex');
    text(0.026,1e-6,'$k^{\rm max}$','FontSize',10,'Interpreter','latex');
    arrow([0.032 3e-6],[0.032 3e-3],5,'BaseAngle',60,'EdgeColor','b','FaceColor','b','LineWidth',1.0);
    %arrow([0.086 1e-2],[0.086 3e-4],5,'BaseAngle',60,'EdgeColor',[232/255 0 0],'FaceColor',[232/255 0 0],'LineWidth',1.0);
    
    subplot(2,1,2);   %0.043851585069891
    
    h(3) = loglog(kx2,zfx2,'b');
    hold on;
    loglog(0.1*[1 1],[1e-8 1e-1],'--k','LineWidth',0.4);
    h(4) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    %loglog(1/sqrt(200)*[1 1],[1e-9 1],'--k','LineWidth',0.5);
    %loglog([1 1],[1e-9 1],'--k','LineWidth',0.5);
    kk = linspace(3.4e-2,9e-2,2);
    loglog(kk,8e-7*kk.^(-3),'k','LineWidth',0.6);
    kk = linspace(1e-1,4e-1,2);
    loglog(kk,1.2e-8*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-7 1e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'YMinorTick','on');
    xlabel('$k$');
    ylabel('$|\delta\!B_{z\!,k~}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.09,0.84,0.43]);
    plotTickLatex2D('xtickdy',-0.00,'xlabeldy',0.012,'ytickdx',-0.005,'ylabeldx',-0.0285,'ytickdy',0.01);
    text(6e-2,1.2e-2,'$\propto\! k^{-3}$','FontSize',10,'Interpreter','latex');
    text(2.1e-1,8e-5,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(3.1e-3,8e-7,'$St = 1$','FontSize',10,'Interpreter','latex');
    

    
   
     
    
    
    
    
    
    
    
    
elseif i==7
%
% distribution function for MIRROR
%

set(0,'DefaultLineLineWidth',2.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',40);
set(0,'DefaultTextFontSize',40);
set(0,'DefaultAxesLineWidth',2.0);

shear = 3e-4;

f = 2;

if (f<10)
    numlab = ['000',num2str(f)];
else if (f<100)
    numlab = ['00',num2str(f)];
else if (f<1000)
    numlab = ['0',num2str(f)];
else
    numlab = num2str(f);
end
end
end
 
fid=fopen(['/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/par-alt/fhshear.',numlab,'.par.lis'],'rb');
%fid=fopen(['/Volumes/My Passport Studio/pegasus/mrs-slow/par/mrshear.',numlab,'.par.lis'],'rb');

% Read the coordinate limits
coorlim = fread(fid,12,'float');
x1l = coorlim(1);     x1u = coorlim(2);
x2l = coorlim(3);     x2u = coorlim(4);
x3l = coorlim(5);     x3u = coorlim(6);
x1dl = coorlim(7);    x1du = coorlim(8);
x2dl = coorlim(9);    x2du = coorlim(10);
x3dl = coorlim(11);   x3du = coorlim(12);
 
% Read number of particle types
npartypes = fread(fid,1,'int32');
for i=1:npartypes
   attr = fread(fid,2,'float');
   mass(i) = attr(1);
   qomc(i) = attr(2);
end
 
% Read the time
time = fread(fid,1,'float')*shear
dt   = fread(fid,1,'float');
 
% Read the particle number
n =  fread(fid,1,'int64')
 
% Read all the particle information
xbins = 100; ybins = 60;
xrange=linspace(-25,25,xbins); yrange=linspace(0,30,ybins);
count = zeros(xbins-1,ybins-1);
mbins = 1000;

for m=1:mbins
    
p = 1; m

for i=((m-1)*n/mbins)+1:(m*n/mbins)
    parinfo = fread(fid,10,'float');
    %v1(p) = parinfo(4);
    %v2(p) = parinfo(5);
    %v3(p) = parinfo(6);
    %vprl(p) = parinfo(9);
    %vprp(p) = parinfo(10);
    data(p,1) = parinfo(9);
    data(p,2) = parinfo(10);
    prop = fread(fid,1,'int32');
    pid = fread(fid,1,'int64');
    cpuid = fread(fid,1,'int64');
    %vsq(p) = v1(p)*v1(p)+v2(p)*v2(p)+v3(p)*v3(p)1;
    p = p+1;
end

count = count + hist2d(data,xrange,yrange);

end

save('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc220.mat', ...
     'xbins','ybins','xrange','yrange','count');

 %{
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc0.mat');
count0=count;
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc220.mat');
count220=count;
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc330.mat');
count330=count;     
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc500.mat');
count500=count;
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc1000.mat');
count1000=count;
load('/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/distfunc3000.mat');
count3000=count;

for i=1:99
    xx(i) = xrange(i)+0.5;
end
for i=1:59
    yy(i) = yrange(i)+0.5;
end
for i=1:59
    fvp0(:,i) = count0(:,i)*yy(i);
    fvp220(:,i) = count220(:,i)*yy(i);
    fvp330(:,i) = count330(:,i)*yy(i);
    fvp500(:,i) = count500(:,i)*yy(i);
    fvp1000(:,i) = count1000(:,i)*yy(i);
    fvp3000(:,i) = count3000(:,i)*yy(i);
end
pcolor(xx,yy,fvp220'-fvp0'); shading interp; axis equal;

fvprl0 = sum(fvp0,2);
fvprl220 = sum(fvp220,2);
fvprl330 = sum(fvp330,2);
fvprl500 = sum(fvp500,2);
fvprl1000 = sum(fvp1000,2);
fvprl3000 = sum(fvp3000,2);
     
     
     
load('/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/distfunc0.mat');
count0=count;
load('/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/distfunc1.mat');
count1=count;
load('/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/distfunc6.mat');
count6=count;

for i=1:99
    xx(i) = xrange(i)+0.5;
end
for i=1:59
    yy(i) = yrange(i)+0.5;
end
for i=1:59
    fvp0(:,i) = count0(:,i)*yy(i);
    fvp1(:,i) = count1(:,i)*yy(i);
    fvp6(:,i) = count6(:,i)*yy(i);
end
pcolor(xx,yy,fvp6'-fvp0'); shading interp; axis equal;

fvprl0 = sum(fvp0,2);
fvprl1 = sum(fvp1,2);
fvprl6 = sum(fvp6,2);
%}


%{

[nn,xout] = hist(vprl,100);
hold on;
plot(xout,nn/max(nn),'c');
%bar(xout,nn,'barwidth',1,'basevalue',1)
set(gca,'YScale','log')
axis([-40 40 1e-3 1 ])
%}

%vprl = abs(vprl);

%{
xrange = linspace(-48,48,20);
yrange = xrange;
zrange = xrange;
data(:,1) = v1(:);
data(:,2) = v2(:);
data(:,3) = v3(:);

for i=1:length(xrange)-1
    data((data(:,1)>xrange(i))&(data(:,1)<=xrange(i+1)),4)=i;
 end

for i=1:length(yrange)-1
    data((data(:,2)>yrange(i))&(data(:,2)<=yrange(i+1)),5)=i;  
end

for i=1:length(zrange)-1
    data((data(:,3)>zrange(i))&(data(:,3)<=zrange(i+1)),6)=i;
end

count=zeros(length(xrange)-1,length(yrange)-1,length(zrange)-1);

data=data(data(:,4)>0,:); % if a data point is out of the x range, throw it away
data=data(data(:,5)>0,:);% if a data point is out of the y range, throw it away
data=data(data(:,6)>0,:);

for i=1:size(data,1)
    count(data(i,4),data(i,5),data(i,6))=count(data(i,4),data(i,5),data(i,6))+1; 
end

count = count + 1;
flnf = count .* log(count);
ent = sum(sum(sum(flnf)))

%}


%{
fig(f+1);clf;

set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,40,40]);
set(gcf,'PaperPositionMode','auto');

imagesc(xrange,yrange,count');
axis xy; axis image; shading flat;
axis([0 50 0 50]);

%}



elseif i==9
%
% particles tracking for MIRROR
%


fn1=5453;    % first figure index
fn2=6280;  % last figure index

Lx = 1152; Ly = 1152; 
shear = 3e-4;
dir ='/Volumes/My Passport Studio/pegasus/mrs-slow/track/';
%shear = 3e-3;
%dir = '/Volumes/My Passport Studio/pegasus/mrs-fast/track/';        % directory
fname = 'mrshear';

for ff=fn1:fn2
    ff
 if (ff<10)
     numlab = ['000',num2str(ff)];
     else if (ff<100)
         numlab = ['00',num2str(ff)];
         else if (ff<1000)
             numlab = ['0',num2str(ff)];
             else
                 numlab = num2str(ff);
             end
         end
 end
 
 fid=fopen([dir,fname,'.',numlab,'.track.lis'],'rb');
  
 % Read the coordinate limits
 coorlim = fread(fid,12,'float');
 x1l = coorlim(1); x1u = coorlim(2); x2l = coorlim(3); x2u = coorlim(4);
 
 % Read number of particle types
 npartypes = fread(fid,1,'int32');
 for i=1:npartypes
     m = fread(fid,1,'float');
     qomc = fread(fid,1,'float');
 end
 
 % Read the time
 t(ff) = fread(fid,1,'float');
 dt   = fread(fid,1,'float');
 
 % Read the particle number
 npar =  fread(fid,1,'int64');
 
 % Read all the particle information
 for j=1:npar
     parinfo = fread(fid,10,'float');
     %x1(j) = parinfo(1); x2(j) = parinfo(2);
     %v1(j) = parinfo(4); v2(j) = parinfo(5); v3(j) = parinfo(6);
     mu(j) = parinfo(8);
     %vprl(j) = parinfo(9); 
     vprp(j) = parinfo(10);
     prop = fread(fid,1,'int32');
     pid = fread(fid,1,'int64');
     cpuid(j) = fread(fid,1,'int64');
 end
 
 [cpu,ii] = sort(cpuid);
 
 %x1s(ff,:) = x1(ii);
 %x2s(ff,:) = x2(ii);
 %v1s(ff,:) = v1(ii);
 %v2s(ff,:) = v2(ii);
 %v3s(ff,:) = v3(ii);
 mus(ff,:) = mu(ii);
 %vprls(ff,:) = vprl(ii);
 vprps(ff,:) = vprp(ii);
 
 %clear x1 x2 cpuid;
 %clear v1 v2 v3 mu vprl vprp;
 clear mu vprp cpuid;
 
 fclose(fid);
 
end

B = vprps.^2 ./ mus;

it=337;  % 1032, 1814, 2272, 127, 6, 74, 257, 258, 268, 271, 278, 289, 304, 327, 328, 330, 394, 408, 556, 620, 861, 899, 916, 1111, 1126, 1621, 1733
ip=2187; % 146, 531, 1059, 1479

%{

figure(99);clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,11.5,5.5]);
set(gcf,'PaperPositionMode','auto');

subplot(2,1,1);
plot(t*shear,smooth(vprls(:,it),10),'b');
hold on;
plot(t*shear,smooth(vprls(:,ip),10),'Color',[232/255 0 0]);
hold off;
axis([0 1.8 -6 6]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'XTickLabel',['']);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$v_{||}$');
set(gca,'Units','normalized','Position',[0.095 0.58 0.875 0.41]);
plotTickLatex2D('ylabeldx',-0.005,'ytickdy',0.01,'ytickdx',0.003);
subplot(2,1,2);
plot(t*shear,smooth(mus(:,it),10),'b');
hold on;
plot(t*shear,smooth(mus(:,ip),10),'Color',[232/255 0 0]);
hold off;
axis([0 1.8 0 45]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'YTick',[10 20 30 40]);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$\mu$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.095 0.14 0.875 0.41]);
plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ylabeldx',-0.02,'ytickdy',0.01,'ytickdx',0.003);
%}

figure(9);clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,11.5,8]);
set(gcf,'PaperPositionMode','auto');

subplot(3,1,1);
semilogy(t*shear,smooth(mus(:,it),10)/mus(1,it),'b');
hold on;
semilogy(t*shear,smooth(mus(:,ip),10)/mus(1,ip),'Color',[232/255 0 0]);
semilogy([0.16 0.16],[2e-2 20],'--k','Linewidth',0.5);
semilogy([1.4 1.4],[2e-2 20],'--k','Linewidth',0.5);
hold off;
axis([0 1.8 2e-2 20]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'YTick',[1e-1 1 10]);
set(gca,'YMinorTick','on','XMinorTick','on','XTickLabel',['']);
ylabel('$\mu / \mu_0$');
set(gca,'Units','normalized','Position',[0.1 0.7 0.875 0.285]);
plotTickLatex2D('ylabeldx',0.01,'ytickdy',0.005,'ytickdx',0.003);
subplot(3,1,2);
plot(t*shear,smooth(B(:,it),10),'b');
hold on;
plot(t*shear,smooth(B(:,ip),10),'Color',[232/255 0 0]);
plot([0.16 0.16],[0.5 4.5],'--k','Linewidth',0.5);
plot([1.4 1.4],[0.5 4.5],'--k','Linewidth',0.5);
text(0.1,2.2,'linear','Rotation',90,'Fontsize',10);
text(0.66,3.85,'secular','Fontsize',10);
text(1.44,3.85,'saturated','Fontsize',10);
hold off;
axis([0 1.8 0.5 4.5]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'YTick',[1 2 3 4]);
set(gca,'YMinorTick','on','XMinorTick','on','XTickLabel',['']);
ylabel('$B/B_0$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.1 0.4 0.875 0.285]);
plotTickLatex2D('ylabeldx',-0.033,'ytickdy',0.005,'ytickdx',0.003);
subplot(3,1,3)
plot(t*shear,smooth(vprls(:,it),10),'b');
hold on;
plot(t*shear,smooth(vprls(:,ip),10),'Color',[232/255 0 0]);
plot([0.16 0.16],[-6 6],'--k','Linewidth',0.5);
plot([1.4 1.4],[-6 6],'--k','Linewidth',0.5);
hold off;
axis([0 1.8 -6 6]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8]);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$v_{||}\;/v_{\!{\rm A0}}$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.1 0.1 0.875 0.285]);
plotTickLatex2D('xtickdy',0.002,'xlabeldy',0.005,'ylabeldx',-0.005,'ytickdy',0.005,'ytickdx',0.003);




elseif i==10
%
% particles tracking for FIREHOSE
%


fn1=1;    % first figure index
fn2=3333;  % last figure index

Lx = 1152; Ly = 1152; shear = 3e-4;
dir ='/Volumes/My Passport Studio/pegasus/fhs-slow/track/';        % directory
fname = 'fhshear';

for ff=fn1:fn2
    ff
 if (ff<10)
     numlab = ['000',num2str(ff)];
     else if (ff<100)
         numlab = ['00',num2str(ff)];
         else if (ff<1000)
             numlab = ['0',num2str(ff)];
             else
                 numlab = num2str(ff);
             end
         end
 end
 
 fid=fopen([dir,fname,'.',numlab,'.track.lis'],'rb');
  
 % Read the coordinate limits
 coorlim = fread(fid,12,'float');
 x1l = coorlim(1); x1u = coorlim(2); x2l = coorlim(3); x2u = coorlim(4);
 
 % Read number of particle types
 npartypes = fread(fid,1,'int32');
 for i=1:npartypes
     m(i) = fread(fid,1,'float');
     qomc(i) = fread(fid,1,'float');
 end
 
 % Read the time
 t(ff) = fread(fid,1,'float');
 dt   = fread(fid,1,'float');
 
 % Read the particle number
 npar =  fread(fid,1,'int64');
 
 % Read all the particle information
 for j=1:npar
     parinfo = fread(fid,10,'float');
     x1(j) = parinfo(1); x2(j) = parinfo(2);
     v1(j) = parinfo(4); v2(j) = parinfo(5); v3(j) = parinfo(6);
     mu(j) = parinfo(8);
     vprl(j) = parinfo(9); 
     vprp(j) = parinfo(10);
     prop = fread(fid,1,'int32');
     pid = fread(fid,1,'int64');
     cpuid(j) = fread(fid,1,'int64');
 end
 
 [cpu,ii] = sort(cpuid);
 
 x1s(ff,:) = x1(ii);
 x2s(ff,:) = x2(ii);
 v1s(ff,:) = v1(ii);
 v2s(ff,:) = v2(ii);
 v3s(ff,:) = v3(ii);
 mus(ff,:) = mu(ii);
 vprls(ff,:) = vprl(ii);
 vprps(ff,:) = vprp(ii);
 
 clear x1 x2;
 clear v1 v2 v3 mu vprl vprp;
 clear cpuid;
 
 fclose(fid);
 
end

%{
subplot(2,1,1);
plot(t*shear,smooth(vprls(:,it),10),'b');
hold on;
plot(t*shear,smooth(vprls(:,ip),10),'Color',[232/255 0 0]);
hold off;
axis([0 1.8 -6 6]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'XTickLabel',['']);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$v_{||}$');
set(gca,'Units','normalized','Position',[0.095 0.58 0.875 0.41]);
plotTickLatex2D('ylabeldx',-0.005,'ytickdy',0.01,'ytickdx',0.003);
subplot(2,1,2);
plot(t*shear,smooth(mus(:,it),10),'b');
hold on;
plot(t*shear,smooth(mus(:,ip),10),'Color',[232/255 0 0]);
hold off;
axis([0 1.8 0 45]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8],'YTick',[10 20 30 40]);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$\mu$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.095 0.14 0.875 0.41]);
plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ylabeldx',-0.02,'ytickdy',0.01,'ytickdx',0.003);
%}

B = vprps.^2 ./ mus;

%43,1126,1354,1976,2193

it = 2193;
t = 1:3333;

%{

figure(99);clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,11.5,8]);
set(gcf,'PaperPositionMode','auto');

subplot(3,1,1);
plot(t*shear,smooth(vprls(:,it),10),'k');
axis([0 1 -25 25]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'XTickLabel',['']);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$v_{||}$');
set(gca,'Units','normalized','Position',[0.11 0.7 0.875 0.285]);
plotTickLatex2D('ylabeldx',0,'ytickdy',0.005,'ytickdx',0.003);
subplot(3,1,2);
semilogy(t*shear,smooth(mus(:,it),10)/mus(1,it),'k');
axis([0 1 1e-1 2]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'YTick',[1e-1 1 10]);
set(gca,'YMinorTick','on','XMinorTick','on','XTickLabel',['']);
ylabel('$\mu / \mu_0$');
set(gca,'Units','normalized','Position',[0.11 0.4 0.875 0.285]);
plotTickLatex2D('ylabeldx',0.0,'ytickdy',0.005,'ytickdx',0.003);
subplot(3,1,3);
plot(t*shear,smooth(B(:,it),10),'k');
axis([0 1 0.5 1.1]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'YTick',[0.6 0.8 1]);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$B$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.11 0.1 0.875 0.285]);
plotTickLatex2D('xtickdy',0.002,'xlabeldy',0.005,'ylabeldx',-0.026,'ytickdy',0.005,'ytickdx',0.003);
%}

figure(10);clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,11.5,5.5]);
set(gcf,'PaperPositionMode','auto');

subplot(2,1,1);
semilogy(t*shear,smooth(mus(:,it),10)/mus(1,it),'k');
hold on;
semilogy([0.08 0.08],[1e-1 2],'--k','Linewidth',0.5);
semilogy([0.32 0.32],[1e-1 2],'--k','Linewidth',0.5);
hold off;
axis([0 1 1e-1 2]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'YTick',[1e-1 1 10]);
set(gca,'YMinorTick','on','XMinorTick','on','XTickLabel',['']);
ylabel('$\mu / \mu_0$');
set(gca,'Units','normalized','Position',[0.11 0.585 0.875 0.41]);
plotTickLatex2D('ylabeldx',0.0,'ytickdy',0.005,'ytickdx',0.003);
subplot(2,1,2);
plot(t*shear,smooth(B(:,it),10),'k');
hold on;
plot([0.08 0.08],[0.5 1.1],'--k','Linewidth',0.5);
plot([0.32 0.32],[0.5 1.1],'--k','Linewidth',0.5);
hold off;
axis([0 1 0.5 1.1]);
set(gca,'TickLength',[0.025 0.025],'FontSize',10);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'YTick',[0.6 0.8 1]);
set(gca,'YMinorTick','on','XMinorTick','on');
ylabel('$B/B_0$'); xlabel('$St$');
set(gca,'Units','normalized','Position',[0.11 0.15 0.875 0.41]);
plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.025,'ylabeldx',-0.0185,'ytickdy',0.005,'ytickdx',0.003);
text(0.05,0.59,'linear','Fontsize',10,'Rotation',90);
text(0.135,0.63,'secular','Fontsize',10);
text(0.37,0.63,'saturated','Fontsize',10);



elseif i==11
%
% nueff plot for firehose & mirror
%
%ffas   check! 0.0494 for secular phase
%fmed
%fslo   check! 0.000152, sec: 0.00586 from 1/mean(b); 0.003851 from width; js=200; je=500; avg is 0.00486
%fss    check! sec: 0.00566 from 1/mean(b); 0.00327 from width; js=400; je=1200; avg is 0.0045

%mfas   check!
%mmed
%mslo
%mss

 
    fig(11); clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,4.5]);
    set(gcf,'PaperPositionMode','auto');
    
    %firehose
    subplot(1,2,1);
    loglog(2*0.01,0.017,'+k','Markersize',9);
    hold on;
    loglog(2*0.01,0.017,'sk','Markersize',5);
    loglog(2*0.01,0.00233,'xk','Markersize',9);
    loglog(2*0.03,0.050,'+b','Markersize',9);
    loglog(2*0.03,0.051,'sb','Markersize',5);
    loglog(2*0.03,0.00486,'xb','Markersize',9);
    loglog(2*0.1,0.12,'+','Markersize',9,'Color',[34/255 139/255 34/255]);
    loglog(2*0.1,0.12,'s','Markersize',5,'Color',[34/255 139/255 34/255]);
    loglog(2*0.1,0.0172,'x','Markersize',9,'Color',[34/255 139/255 34/255]);
    loglog(2*0.3,0.26,'+','Markersize',9,'Color',[232/255 0 0]);
    loglog(2*0.3,0.26,'s','Markersize',5,'Color',[232/255 0 0]);
    loglog(2*0.3,0.0494,'x','Markersize',9,'Color',[232/255 0 0]);
    shbeta = linspace(1.5e-2,6.7e-1,2);
    loglog(shbeta,shbeta*0.85,'--k','LineWidth',0.4);
    loglog(1.4e-2,2.6e-1,'+k','Markersize',6);
    text(1.7e-2,2.21e-1,'$\nu_{\rm scatt,sat}$','Fontsize',10);
    loglog(1.7e-1,2e-3,'xk','Markersize',7);
    text(2.1e-1,1.75e-3,'$\nu_{\rm scatt,sec}$','Fontsize',10);
    loglog(1.4e-2,1.4e-1,'sk','Markersize',5);
    text(1.7e-2,1.3e-1,'$\nu_{\rm f}$','Fontsize',10);
    hold off;
    axis([1e-2 1 1e-3 1e0]);
    set(gca,'TickLength',[0.035 0.04],'Fontsize',10);
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel('$S\!\beta_0$'); set(gca,'XTick',[1e-2 1e-1 1e0],'YTick',[1e-2 1e-1 1e0]);
    set(gca,'Units','normalized','Position',[0.075,0.17,0.444,0.79]);
    plotTickLatex2D('ytickdx',-0.00,'ylabeldx',-0.02,'ytickdy',0.01,'xtickdy',-0.02,'xlabeldy',-0.03);
    text(1.3e-2,5e-1,'firehose','Fontsize',10);
    text(0.105,0.33,'$\propto\!S\!\beta_0$','Fontsize',10);

    %mirror
    subplot(1,2,2);
    loglog(2*0.01,0.0091,'+k','Markersize',9);
    hold on;
    loglog(2*0.01,0.0091,'sk','Markersize',5);
    %loglog(2*0.01,0.0039,'+k');
    %loglog(2*0.01,0.0127,'xk');
    loglog(2*0.01,0.0022,'xk','Markersize',9);
    
    loglog(2*0.03,0.028,'+b','Markersize',9);
    loglog(2*0.03,0.029,'sb','Markersize',5);
    %loglog(2*0.03,0.004,'+b');
    %loglog(2*0.03,0.013,'xb');
    loglog(2*0.03,0.0021,'xb','Markersize',9);
    
    loglog(2*0.1,0.092,'+','Markersize',9,'Color',[34/255 139/255 34/255]);
    loglog(2*0.1,0.094,'s','Markersize',5,'Color',[34/255 139/255 34/255]);
    %loglog(2*0.1,0.004,'+','Color',[34/255 139/255 34/255]);
    %loglog(2*0.1,0.0128,'x','Color',[34/255 139/255 34/255]);
    loglog(2*0.1,0.0022,'x','Markersize',9,'Color',[34/255 139/255 34/255]);
    
    loglog(2*0.3,0.055,'+','Markersize',9,'Color',[232/255 0 0]);
    loglog(2*0.3,0.054,'s','Markersize',5,'Color',[232/255 0 0]);
	loglog(2*0.3,0.011,'x','Markersize',9,'Color',[232/255 0 0]);


    shbeta = linspace(1.5e-2,6.7e-1,2);
    loglog(shbeta,shbeta*0.46,'--k','LineWidth',0.4);
    loglog(1.4e-2,2.6e-1,'+k','Markersize',6);
    text(1.7e-2,2.21e-1,'$\nu_{\rm scatt,sat}$','Fontsize',10);
    loglog(1.4e-2,1.47e-1,'xk','Markersize',7);
    text(1.7e-2,1.3e-1,'$\nu_{\rm scatt,sec}$','Fontsize',10);
    loglog(1.4e-2,0.83e-1,'sk','Markersize',5);
    text(1.7e-2,0.767e-1,'$\nu_{\rm m}$','Fontsize',10);
    hold off;
    axis([1e-2 1 1e-3 1e0]);
    set(gca,'XMinorTick','on','YMinorTick','on','YTickLabel',['']);
    set(gca,'TickLength',[0.035 0.04],'Fontsize',10);
    set(gca,'Units','normalized','Position',[0.530,0.17,0.444,0.79]);
    xlabel('$S\!\beta_0$'); set(gca,'XTick',[1e-1 1e0],'YTick',[1e-3 1e-2 1e-1 1e0])
    plotTickLatex2D('xtickdy',-0.02,'xlabeldy',-0.03);
    text(1.28e-2,5e-1,'mirror','Fontsize',10);
    text(0.12,0.21,'$\propto\!S\!\beta_0$','Fontsize',10);
    text(0.036,0.006,'trapped','Fontsize',8,'Color','k');
    drawbrace([1.66e-2,2.5e-3],[2.4e-1,2.5e-3],.12,'Color','k','Linewidth',0.4);
    %36% trapped in super saturated
    
%sup: secular (St = 0.12), 2e-5, 0.0127 (0.0022 trapped, 0.038 passing), 71% trapped; saturated (St = ), 3.7e-5, 0.0091 ( 0.0039 trapped, 0.0244 passing)
%slow: secular (St = 0.2 - 1.2), 7.9e-5, 0.013 (0.0021 trapped, 0.033 passing); saturated (St = 1.66 - 1.9), 0.00017, 0.028 (0.004 trapped, 0.04 passing)
%med: 
%fast: secular (St = ), 48.7e-5, ( trapped, passing); saturated (St = ), 97.4e-5, 0.052 (trapped, passing)


%{

    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,5.5]);
    set(gcf,'PaperPositionMode','auto');
    

    loglog(0.02,-0.009,'+k');
    hold on;
    text(0.0187,0.0087,'${\rm m}_{\rm sat}$','Fontsize',10,'Color','k');
    %loglog(0.02,0.0039,'+k');
    %loglog(0.02,0.0127,'xk');
    %loglog(0.02,0.0022,'xk');
    text(0.0187,0.0021,'${\rm m}^{\rm trap}_{\rm sec}$','Fontsize',10,'Color','k');
    %loglog(0.06,0.028,'+b');
    text(0.0562,0.0265,'${\rm m}_{\rm sat}$','Fontsize',10,'Color','b');
    %loglog(0.06,0.004,'+b');
    %loglog(0.06,0.013,'xb');
    %loglog(0.06,0.0021,'xb');
    text(0.0562,0.00201,'${\rm m}^{\rm trap}_{\rm sec}$','Fontsize',10,'Color','b');
    %loglog(0.2,0.09,'+','Color',[34/255 139/255 34/255]);
    text(0.187,0.086,'${\rm m}_{\rm sat}$','Fontsize',10,'Color',[34/255 139/255 34/255]);
    %loglog(0.2,0.004,'+','Color',[34/255 139/255 34/255]);
    %loglog(0.2,0.0128,'x','Color',[34/255 139/255 34/255]);
    %loglog(0.2,0.0022,'x','Color',[34/255 139/255 34/255]);
    text(0.187,0.0022,'${\rm m}^{\rm trap}_{\rm sec}$','Fontsize',10,'Color',[34/255 139/255 34/255]);
    %loglog(0.6,0.34,'+','Color',[232/255 0 0]);
    text(0.562,0.324,'${\rm m}_{\rm sat}$','Fontsize',10,'Color',[232/255 0 0]);
    
    text(0.08,0.19,'$\propto\!S \beta$','Fontsize',10);
    shbeta = linspace(1.5e-2,6.7e-1,2);
    loglog(shbeta,shbeta*0.85,'--k','LineWidth',0.4);
    shbeta = linspace(1.5e-2,6.7e-1,2);
    loglog(shbeta,shbeta*0.46,'--k','LineWidth',0.4);
    
    %loglog(0.02,0.017,'+k');
    text(0.01956,0.017,'${\rm f}_{\rm sat}$','Fontsize',10,'Color','k');
    %loglog(0.06,0.051,'+b');
    text(0.0587,0.051,'${\rm f}_{\rm sat}$','Fontsize',10,'Color','b');
    %loglog(0.2,0.16,'+','Color',[34/255 139/255 34/255]);
    text(0.1956,0.16,'${\rm f}_{\rm sat}$','Fontsize',10,'Color',[34/255 139/255 34/255]);
    %loglog(0.6,0.34,'+','Color',[232/255 0 0]); 
    text(0.587,0.34,'${\rm f}_{\rm sat}$','Fontsize',10,'Color',[232/255 0 0]);
        
    hold off;
    axis([1e-2 1 1e-3 1e0]);
    set(gca,'TickLength',[0.025 0.025],'Fontsize',10);
    ylabel('$\nu_{\rm e\!ff}$'); xlabel('$S \beta$');
    set(gca,'XMinorTick','on','YMinorTick','on');
    %set(gca,'YTick',[-0.1 0 0.1 0.2 0.3]); 
    %set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gca,'Units','normalized','Position',[0.13,0.15,0.845,0.815]);
    plotTickLatex2D('xtickdy',-0.015,'xlabeldy',-0.02,'ytickdx',-0.0,'ylabeldx',-0.02,'ytickdy',0.01);

    %}
    
elseif i==55
    
    
Lx = 1152; Ly = 1152; shear = 3e-4;

%dir ='/Volumes/My Passport Studio/pegasus/fhs-super-slow/vtk/';
dir ='/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/vtk/';        % directory
fname = 'fhshear';

    f=333;
        
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
    end

    % declare file name
    filename = [dir,fname,'.',numlab,'.vtk'];

    % open file and initialize grid
    ary = [3 1 3 9];
    [Grid,status] = init_grid(filename,ary);

    [time,var,name,status] = readvtk(Grid,filename,1);
    U.bz  = squeeze(var(3,:,:,:));
    
    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d   = squeeze(var);
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = 3/sqrt(13);
    bx  = 2/sqrt(13);
    by  = by0 - bx*shear*time;
    dBz = U.bz;
           
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    end

    zfx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
    end
    zfx1 = zfx/nsy;
    %scale = sum(sum((abs(fft2(dBz))).^2))/(Grid.nx1*Grid.nx2)^2/sum(zfx1);
    %zfx1 = zfx1*scale;
    
    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx1 = kx(1:nsx/2) * (2*pi/pstep);

% Calculate the frequencies

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx2;			% step taken in y direction during field-line trace
    nsy = 2*Ly/Grid.dx2;			% number of steps taken to get periodic series
    nsx = Lx/Grid.dx1;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBz',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    end
    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5);
    zfy1 = sum(zfy,2)/nsx/nsy;
    ky1 = ky*2*pi;
    
    %{
    zfy = zeros(nsy/2,1);
    for n=1:nsx
        zzf = fft(zz(:,n))/nsy;
        zfy = zfy + 2*abs(zzf(1:nsy/2)).^2;
    end
    zfy1 = zfy/nsx;
    scale = sum(sum((abs(fft2(dBz))).^2))/(Grid.nx1*Grid.nx2)^2/sum(zfy1);
    zfy1 = zfy1*scale;
    %}
    %{
    ky = mod( 1/2 + (0:(nsy-1))/nsy , 1 ) - 1/2;
    ky1 = ky(1:nsy/2) * (2*pi/pstep);
    %}
  

    img = U.d-mean(mean(U.d));    
    [N M] = size(img);
    imgf = fftshift(fft2(img));
    imgfp = (abs(imgf)/(N*M)).^2;
    dimDiff = abs(N-M);
    dimMax = max(N,M);
    % Make square
    if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
    elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
    end

    halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)
    [X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
    [theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
    rho = round(rho);
    i = cell(floor(dimMax/2) + 1, 1);
    for r = 0:floor(dimMax/2)
        i{r + 1} = find(rho == r);
    end
    Pf = zeros(1, floor(dimMax/2)+1);
    for r = 0:floor(dimMax/2)
        Pf(1, r + 1) = nansum( imgfp( i{r+1} ) );
    end
    f1 = linspace(0,halfDim-1,length(Pf));                                           % Set abscissa
    ky3 = f1*2*pi/1152.;
    zfy3 = Pf./ky3.^2;%/(2*pi*ky3).^2;
    
    
    
    
    
Lx = 1152; Ly = 1152; shear = 3e-4;
dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/vtk/';        % directory
fname = 'mrshear';
    
   f = 333;
    
    % format file number
    if (f<10)
        numlab = ['000',num2str(f)];   
    elseif (f<100)
       numlab = ['00',num2str(f)];
    elseif (f<1000)
       numlab = ['0',num2str(f)];
    else
       numlab = num2str(f);
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
    
    [time,var,name,status] = readvtk(Grid,filename,2);
    U.d   = squeeze(var); dn = U.d - mean(mean(U.d));
    
    for j=1:Grid.nx1
        yy(j,:) = Grid.x2(:);
    end
    for i=1:Grid.nx2
        xx(:,i) = Grid.x1(:);
    end
    by0 = -1/sqrt(5);
    bx  = 2/sqrt(5);
    by  = by0 - bx*shear*time;
    dBy = U.by-by;
    dBx = U.bx-bx;
    dBprl = (by*dBy+bx*dBx)/sqrt(by*by+bx*bx);
    Bmag = sqrt(U.bx.^2+U.by.^2+U.bz.^2);
    
    %parallel spectrum
    a = by0/bx - shear*time;		% instantaneous pitch angle of mean magnetic field
    x0 = -0.5*Lx;                   % where to start in x
    y0 = -0.5*Ly;                   % where to start in y
    x0step =  0;                    % step taken in x direction for x0
    y0step =  Grid.dx2;             % step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);	% step taken in parallel direction during field-line trace
    xstep = Grid.dx1;               % step taken in x direction during field-line trace
    ystep = a*Grid.dx1;             % step taken in y direction during field-line trace
    nsx = 2*Lx/Grid.dx1;			% number of steps taken to get periodic series
    nsy =   Ly/Grid.dx2;			% number of steps taken to shift y0 and cover space
    zz = zeros(nsx,nsy); %zn = zeros(nsx,nsy);

    for n=1:nsy
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsx+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,Bmag',xi,yi);
    %zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsx,n) = inpaint_nans(zi(1:nsx));

    %ni = interp2(Grid.x1,Grid.x2,dn',xi,yi);
    %zn(1:nsx,n) = inpaint_nans(ni(1:nsx));
    end

    zfx = zeros(nsx/2,1); %znx = zeros(nsx/2,1);
    for n=1:nsy
        zzf = fft(zz(:,n))/nsx;
        zfx = zfx + 2*abs(zzf(1:nsx/2)).^2;
        %znf = fft(zn(:,n))/nsx;
        %znx = znx + 2*abs(znf(1:nsx/2)).^2;
    end
    zfx2 = zfx/nsy;

    kx = mod( 1/2 + (0:(nsx-1))/nsx , 1 ) - 1/2;
    kx2 = kx(1:nsx/2) * (2*pi/pstep);

    %perpendicular spectrum
    x0 =  0.5*Lx;				% where to start in x
    y0 = -0.5*Ly;				% where to start in y
    x0step = -Grid.dx1;			% step taken in x direction for x0
    y0step =  0;				% step taken in y direction for y0
    pstep =  Grid.dx1*sqrt(1+a^2);		% step taken in perpendicular direction for (x0,y0)
    xstep = -a*Grid.dx1;			% step taken in x direction during field-line trace
    ystep =    Grid.dx1;			% step taken in y direction during field-line trace
    nsy =  2*Ly/Grid.dx1;			% number of steps taken to get periodic series
    nsx =  Lx/Grid.dx1;			% number of steps taken to shift x0 and cover space
    zz = zeros(nsy,nsx);% zn = zeros(nsy,nsx);

    for n=1:nsx
    clear xi yi;

    x0 = x0 + x0step;
    y0 = y0 + y0step;
    if (x0 > 0.5*Lx)
        x0 = x0-Lx;
        y0 = y0+shear*Lx*time;
    end
    if (x0 < -0.5*Lx)
        x0 = x0+Lx;
        y0 = y0-shear*Lx*time;
    end
    if (y0 > 0.5*Ly)
        y0 = y0-Ly;
    end
    if (y0 < -0.5*Ly)
        y0 = y0+Ly;
    end

    xi(1) = x0;				% update starting location
    yi(1) = y0;

    for i=2:nsy+2			% trace field line in shearing-periodic box
        xi(i) = xi(i-1)+xstep;
        yi(i) = yi(i-1)+ystep;
        if (xi(i) > 0.5*Lx)
            xi(i) = xi(i)-Lx;
            yi(i) = yi(i)+shear*Lx*time;
        end
        if (xi(i) < -0.5*Lx)
            xi(i) = xi(i)+Lx;
            yi(i) = yi(i)-shear*Lx*time;
        end
        if (yi(i) > 0.5*Ly)
            yi(i) = yi(i)-Ly;
        end
        if (yi(i) < -0.5*Ly)
            yi(i) = yi(i)+Ly;
        end
    end

    zi = interp2(Grid.x1,Grid.x2,dBprl',xi,yi);
    zz(1:nsy,n) = inpaint_nans(zi(1:nsy));

    %ni = interp2(Grid.x1,Grid.x2,dn',xi,yi);
    %zn(1:nsy,n) = inpaint_nans(ni(1:nsy));
    end

    
    [zfy,ky,Conf] = psdtsh(zz,pstep,2,1.5); %[zny,ky,Conf] = psdtsh(zn,pstep,2,1.5);
    zfy2 = sum(zfy,2)/nsy/nsx;% zny2 = sum(zny,2)/nsy/nsx;
    ky2 = ky*2*pi;
    
    %{
    zfy = zeros(nsy/2,1);
    for n=1:nsx
        zzf = fft(zz(:,n))/nsy;
        zfy = zfy + 2*abs(zzf(1:nsy/2)).^2;
    end
    zfy2 = zfy/nsx;
    ky = mod( 1/2 + (0:(nsy-1))/nsy , 1 ) - 1/2;
    ky2 = ky(1:nsy/2) * (2*pi/pstep);
    %}
    
    img = dn;
    [N M] = size(img);
    imgf = fftshift(fft2(img));
    imgfp = (abs(imgf)/(N*M)).^2;
    dimDiff = abs(N-M);
    dimMax = max(N,M);
    % Make square
    if N > M                                                                % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
    elseif N < M                                                            % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
    end

    halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)
    [X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
    [theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
    rho = round(rho);
    i = cell(floor(dimMax/2) + 1, 1);
    for r = 0:floor(dimMax/2)
        i{r + 1} = find(rho == r);
    end
    Pf = zeros(1, floor(dimMax/2)+1);
    for r = 0:floor(dimMax/2)
        Pf(1, r + 1) = nansum( imgfp( i{r+1} ) );
    end
    f1 = linspace(0,halfDim-1,length(Pf));                                           % Set abscissa
    kx3 = f1*2*pi/1152.;
    zfx3 = Pf./kx3.^2;%/(2*pi*kx3).^2;
    

    figure(55);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,8.2]);
    set(gcf,'PaperPositionMode','auto');
    
    subplot(2,1,1);

    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.1*[1 1],[1e-8 1e-1],'--k','LineWidth',0.4);
    kk = linspace(3.4e-2,9e-2,2);
    loglog(kk,8e-7*kk.^(-3),'k','LineWidth',0.6);
    kk = linspace(1e-1,4e-1,2);
    loglog(kk,1.2e-8*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-7 1e-1]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on','XTickLabel',['']);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'YMinorTick','on');
    ylabel('$|\delta\!B_{z\,}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.54,0.84,0.43]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.031,'ytickdy',0.01);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$~k_{||}$','$~k_\perp$','$~\rho_{\rm i}^{-1}$','Location','Southwest');
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(6e-2,1.2e-2,'$\propto\! k^{-3}$','FontSize',10,'Interpreter','latex');
    text(2.1e-1,8e-5,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(0.56,1.5e-2,'firehose','FontSize',10);
    text(0.57,3.3e-3,'$St = 1$','FontSize',10,'Interpreter','latex');
    text(0.33,1.3e-2,'(a)','Fontsize',10);
    
    subplot(2,1,2);
    
    h(4) = loglog(kx2,zfx2,'b');
    hold on;
    h(5) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    h(6) = loglog(0.107*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    kk = linspace(3.5e-2,1.07e-1,2);
    loglog(kk,3.2e-8*kk.^(-11/3),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,2.5e-9*kk.^(-4.8),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,4e-1,2);
    loglog(kk,2.2e-7*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 1e-6 1e0]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0],'YMinorTick','on');
    xlabel('$k$');
    ylabel('$|\delta\!B_{||\,}/B_0|^2$');
    set(gca,'Units','normalized','Position',[0.15,0.095,0.84,0.43]);
    plotTickLatex2D('xtickdy',-0.0,'xlabeldy',0.006,'ytickdx',-0.005,'ylabeldx',-0.028,'ytickdy',0.01);
    text(1.8e-2,3e-4,'$\propto\! k^{-11/3}_{||}$','FontSize',10,'Interpreter','latex');
    text(4e-2,9e-6,'$\propto\! k^{-4.8}_{||}$','FontSize',10,'Interpreter','latex');
    text(2.1e-1,1.6e-3,'$\propto\! k^{-4.8}_\perp$','FontSize',10,'Interpreter','latex');
    text(0.56,1.4e-1,'mirror','FontSize',10);
    text(0.57,3.3e-2,'$St = 1$','FontSize',10,'Interpreter','latex');
    text(0.33,1.3e-1,'(b)','Fontsize',10);
    
    figure(555);clf;
    set(gcf,'Color',[1 1 1]);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[1,1,11.5,11.85]);
    set(gcf,'PaperPositionMode','auto');
    
    subplot(3,1,1);

    h(1) = loglog(kx1,zfx1,'b');
    hold on;
    h(2) = loglog(ky1,zfy1,'Color',[232/255 0 0]); 
    h(3) = loglog(0.1*[1 1],[5e-8 5e-2],'--k','LineWidth',0.4);
    kk = linspace(3.2e-2,9.5e-2,2);
    loglog(kk,3.5e-7*kk.^(-3),'k','LineWidth',0.6);
    kk = linspace(1e-1,4e-1,2);
    loglog(kk,0.55e-8*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 5e-8 5e-2]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on','XTickLabel',['']);
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2],'YMinorTick','on');
    ylabel('$|\delta\!B_{z\,}|^2/B_0^2$');
    set(gca,'Units','normalized','Position',[0.15,0.687,0.84,0.296]);
    plotTickLatex2D('ytickdx',-0.005,'ylabeldx',-0.031,'ytickdy',0.01);
    ah1 = gca;
    leg1 = legend(ah1,h(1:3),'$\,k_{||}$','$\,k_\perp$','$\,\rho_{\rm i}^{-1}$','Location','Southwest');
    legend(ah1,'boxoff');
    set(leg1,'interpreter','latex');
    text(4.7e-2,1.05e-2,'$\propto\! k^{-3}$','FontSize',10,'Interpreter','latex');
    text(2.05e-1,4e-5,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(0.9*0.56,0.75e-2,'firehose','FontSize',10);
    text(0.9*0.57,1.65e-3,'$St = 1$','FontSize',10,'Interpreter','latex');
    text(0.9*0.33,0.65e-2,'(a)','Fontsize',10);
    
    subplot(3,1,2);
    
    h(4) = loglog(kx2,zfx2,'b');
    hold on;
    h(5) = loglog(ky2,zfy2,'Color',[232/255 0 0]); 
    h(6) = loglog(0.107*[1 1],[1e-7 1],'--k','LineWidth',0.4);
    kk = linspace(3e-2,1.07e-1,2);
    loglog(kk,1.5e-8*kk.^(-11/3),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,1.15e-9*kk.^(-4.8),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,4e-1,2);
    loglog(kk,1.05e-7*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 5e-7 0.5]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on','XTickLabel',['']);
    set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'YMinorTick','on');
    ylabel('$|\delta\!B_{||\,}|^2/B_0^2$');
    ah2 = gca;
    leg2 = legend(ah2,h(4:6),'$\,k_{||}$','$\,k_\perp$','$\,\rho_{\rm i}^{-1}$','Location','Southwest');
    legend(ah2,'boxoff');
    set(gca,'Units','normalized','Position',[0.15,0.377,0.84,0.297]);
    plotTickLatex2D('xtickdy',-0.0,'xlabeldy',0.006,'ytickdx',-0.005,'ylabeldx',-0.028,'ytickdy',0.01);
    text(1.8e-2,1.5e-4,'$\propto\! k^{-11/3}_{||}$','FontSize',10,'Interpreter','latex');
    text(4.1e-2,4.5e-6,'$\propto\! k^{-4.8}_{||}$','FontSize',10,'Interpreter','latex');
    text(2.05e-1,0.8e-3,'$\propto\! k^{-4.8}_\perp$','FontSize',10,'Interpreter','latex');
    text(0.9*0.56,0.7e-1,'mirror','FontSize',10);
    text(0.9*0.57,1.65e-2,'$St = 1$','FontSize',10,'Interpreter','latex');
    text(0.9*0.33,0.65e-1,'(b)','Fontsize',10);
    
    subplot(3,1,3);
    
    h(7) = loglog(ky3,zfy3,'Color',[1 128/255 40/255]);
    hold on;
    h(8) = loglog(0.1*[1 1],[2e-8 2e-2],'--','Color',[1 128/255 40/255],'LineWidth',0.4);
    h(9) = loglog(kx3,zfx3,'Color',[34/255 139/255 34/255]);
    h(10)= loglog(0.107*[1 1],[2e-8 2e-2],'--','Color',[34/255 139/255 34/255],'LineWidth',0.4);
    kk = linspace(3e-2,9.5e-2,2);
    loglog(kk,40*2.4e-11*kk.^(-3),'k','LineWidth',0.6);
    kk = linspace(1e-1,2.5e-1,2);
    loglog(kk,40*3.9e-13*kk.^(-4.8),'k','LineWidth',0.6);
    kk = linspace(2.7e-2,1.07e-1,2);
    loglog(kk,40*1.1e-10*kk.^(-11/3),'k','LineWidth',0.6);
    kk = linspace(1.1e-1,2.5e-1,2);
    loglog(kk,40*8.7e-12*kk.^(-4.8),'k','LineWidth',0.6);
    hold off;
    axis([2e-3 2e0 2e-8 0.02]);
    set(gca,'TickLength',[0.025 0.025],'FontSize',10);
    set(gca,'XTick',[1e-2 1e-1 1e0],'XMinorTick','on');
    set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2],'YMinorTick','on');
    xlabel('$k d_{\!\rm i0}$');
    ylabel('$|\delta\!n_{\!\rm i}|^2\, k^{-\!2}$');
    ah2 = gca;
    leg2 = legend(ah2,h(7:10),'$\,k\!$,\,firehose','$\,\rho_{\rm i}^{-1}\!$,\,firehose','$\,k\!$,\,mirror','$\,\rho_{\rm i}^{-1}\!$,\,mirror','Location','Southwest');
    legend(ah2,'boxoff');
    set(gca,'Units','normalized','Position',[0.15,0.0675,0.84,0.296]);
    plotTickLatex2D('xtickdy',0.008,'xlabeldy',0.028,'ytickdx',-0.005,'ylabeldx',-0.031,'ytickdy',0.006);
    text(2.3e-2,40*3e-6,'$\propto\! k^{-3}$','FontSize',10,'Interpreter','latex');
    text(4e-2,40*3e-9,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(3.4e-2,40*7e-5,'$\propto\! k^{-11/3}$','FontSize',10,'Interpreter','latex');
    text(1.6e-1,40*2e-7,'$\propto\! k^{-4.8}$','FontSize',10,'Interpreter','latex');
    text(0.9*0.56,40*7.9e-5,'$k$-shell-','FontSize',10);
    text(0.9*0.56,40*1.92e-5,'averaged','FontSize',10);
    text(0.9*0.57,40*4.80e-6,'$St = 1$','FontSize',10,'Interpreter','latex');
    text(0.9*0.33,40*7e-5,'(c)','Fontsize',10);
    

        
else
    error('no plot chosen');
end
