clear all;
figlab = 1;
figM = figure(figlab); hold off;

b = 300; e = 500;

set(0,'DefaultLineLineWidth',1.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultTextFontSize',11);
set(0,'DefaultAxesLineWidth',1.0);

dir ='/Users/kunz/Documents/codes/pegasus/bin/firehose-angle-50/';        % directory
fname = 'firehose';
filename = [dir,fname,'.hst'];

fid = fopen(filename);
C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);

time = C{1};
ME1 = C{22};
ME2 = C{23};
ME3 = C{13};



set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,10,8]);
set(gcf,'PaperPositionMode','auto');

plot(time(b:e),log(ME1(b:e)),'b');
hold on;
plot(time(b:e),log(ME2(b:e)),'r');
plot(time(b:e),log(ME3(b:e)),'g');
%semilogy(time,ME2+ME3,'k');
%semilogy(time,8e-13*exp(2*1.1e-02*time),'--k');
hold off;
axis([0 1000 -35 -2]);


xlabel('$t$');
ylabel('magnetic energy');
set(gca,'TickLength',[0.03 0.03]);
%set(gca,'XTick',[1e1 1e2 1e3]);
%set(gca,'YTick',[1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'Units','normalized','Position',[0.18,0.15,0.75,0.7]);

%{
leg3 = legend('$\delta B_y^2$','$\delta B_z^2$','$\delta B^2$','Location',[0.67    0.2    0.229    0.2333]);
legend(leg3,'boxoff');
set(leg3,'interpreter','latex');
%}

plotTickLatex2D('xtickdy',-0.01,'xlabeldy',-0.02,'ylabeldx',-0.015,'ytickdx',-0.01);
