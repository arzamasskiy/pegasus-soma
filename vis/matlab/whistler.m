set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);

hold off;

ka = logspace(-2,2,128);

for ii=1:128
  waR(ii) = 0.5*ka(ii)^2*( 1.0+sqrt(1.0+(2.0/ka(ii))^2));
  waL(ii) = 0.5*ka(ii)^2*(-1.0+sqrt(1.0+(2.0/ka(ii))^2));
end

figure(1); clf;
subplot(1,2,1);
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,16.5,7]);
set(gcf,'PaperPositionMode','auto');

%for k = [1/16 1/8 ... 8 16], 64 ppc, 64 cells per wavelength
km1  = [1/16   1/8    1/4    1/2    1      2       4       8       16];
wmR1 = [0.0645 0.1331 0.2832 0.6404 1.6180 4.82843 16.9443 64.9848 256.996];
wmL1 = [0.0605 0.1174 0.2207 0.3904 0.6180 0.82843 0.94427 0.98485 0.99612];

%for k = [1 2 4 8 16], 64 ppc, 64 cells in L = 2pi
km2  = [1       2       4       8       16];
wmR2 = [1.61803 4.82843 16.6939 61.3670 205.689];
wmL2 = [0.61803 0.82843 0.92042 0.86159 0.74103];

h = loglog(ka,waR,'-k');
hold on;
loglog(ka,waL,'-k');
f = loglog(km1,wmR1,'+r');
set(f,'MarkerSize',8);
fL = loglog(km1,wmL1,'+r');
set(fL,'MarkerSize',8);
g = loglog(km2,wmR2,'xb');
set(g,'MarkerSize',10);
gL = loglog(km2,wmL2,'xb');
set(gL,'MarkerSize',10);
hold off;
axis([0.02 40 0.02 1e3]);
%axis equal;

xlabel('$k$');
ylabel('$\omega$');
set(gca,'TickLength',[0.03 0.03]);
set(gca,'XTick',[1e-1 1e0 1e1]);
set(gca,'YTick',[1e-1 1e0 1e1 1e2 1e3]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'Units','normalized','Position',[0.08,0.15,0.38,0.82]);

leg1 = legend([f,g,h],'$\,k\Delta x = \pi/32$','$\,\Delta x = \pi/32$','\,theory', ...
    'Location',[0.0887    0.74    0.2245    0.10]);
legend(leg1,'boxoff');
set(leg1,'interpreter','latex');

plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ytickdx',-0.0005,'ylabeldx',-0.01);




subplot(1,2,2);
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,16.5,6]);
set(gcf,'PaperPositionMode','auto');


%for k = 1, 64 ppc, N#d cells per wavelength
N1d   = [4 8 16 32 64 128 256 512 1024];
err1d = [2.075267e-06 8.390385e-07 2.331182e-07 5.062892e-08 1.260410e-08 ...
    3.435530e-09 8.451529e-10 1.970724e-10 6.119160e-11];

N2d = [4 8 16 32 64 128 256];
err2d = [1.833524e-06 9.116288e-07 2.560719e-07 6.109934e-08 1.454531e-08 3.774291e-09 8.603411e-10];

N3d = [4 8 16 32 64 128];
err3d = [1.44325e-06 1.12414e-06 3.44818e-07 7.81062e-08 1.61429e-08 4.61845e-09];

f = loglog(N1d,err1d,'-+k');
hold on;
%f = loglog(N1d,err1d,'+k');
set(f,'MarkerSize',8);

g = loglog(N2d,err2d,'-+b');
set(g,'MarkerSize',8);

h = loglog(N3d,err3d,'-+r');
set(h,'MarkerSize',8);

Na = logspace(1.1,3.1,2);
loglog(Na,0.0004*Na.^(-2),'--k');
hold off;
axis([2 2e3 2e-11 1e-5]);
%axis equal;

xlabel('$N$');
ylabel('$L1$ error');
set(gca,'TickLength',[0.03 0.03]);
set(gca,'XTick',[1e1 1e2 1e3]);
set(gca,'YTick',[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'Units','normalized','Position',[0.6,0.15,0.38,0.82]);

text(100,1.5e-7,'$N^{-2}$');

leg3 = legend('1D','2D','3D','Location',[0.6122    0.1677    0.1304    0.2412]);
legend(leg3,'boxoff');
set(leg3,'interpreter','latex');


plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ylabeldx',-0.01,'ytickdx',-0.001);
