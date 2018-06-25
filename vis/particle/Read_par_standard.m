clear all; clear figM;




set(0,'DefaultLineLineWidth',2.0);
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',40);
set(0,'DefaultTextFontSize',40);
set(0,'DefaultAxesLineWidth',2.0);







bname='mrshear'; % base name
ename='par'; % end name
fn1=11;    % first figure index
fn2=11;  % last figure index

%hold on;
hold off;

tracked = -1;

for f=fn1:fn2
    
    ff = f-fn1+1;

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
 
    %fid=fopen(['/Users/kunz/Documents/codes/pegasus/bin/comb_lis/',bname,'.',numlab,'.',ename,'.lis'],'rb');
    fid =fopen(['/Volumes/My Passport Studio/pegasus/mrs-slow/par/',bname,'.',numlab,'.',ename,'.lis'],'rb');
    %fid=fopen(['/Users/kunz/Documents/codes/pegasus/bin/fhs-slow/par-alt/',bname,'.',numlab,'.',ename,'.lis'],'rb');
    
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
    time(ff) = fread(fid,1,'float')
    dt      = fread(fid,1,'float');
 
    % Read the particle number
    n =  fread(fid,1,'int64')
 
    % Read all the particle information
    p = 1;
    for i=1:100:n
        parinfo = fread(fid,10,'float');
        %if (parinfo(1) > 1200)
         %x1(p) = parinfo(1);
         %x2(p) = parinfo(2);
         %x3(p) = parinfo(3);
         %v1(p) = parinfo(4);
         %v2(p) = parinfo(5);
         %v3(p) = parinfo(6);
         %f0(p) = parinfo(7);
         mu11(p) = parinfo(8);
         %vprl(p) = parinfo(9);
         %vprp(p) = parinfo(10);
         prop = fread(fid,1,'int32');
         pid = fread(fid,1,'int64');
         cpuid(p) = fread(fid,1,'int64');
         %vsq(p) = v1(p)*v1(p)+v2(p)*v2(p)+v3(p)*v3(p);
         p = p+1;
        %end
        
        %{
        if (cpuid(i)==tracked)
            plot(x1(i),x2(i),'.');
            axis([x1l x1u x2l x2u]);
            clear x1 x2 x3 v1 v2 v3 f0 mu vprl vprp prop pid cpuid vmag;
        end
        %}
    end
    
   %[cpu,ii] = sort(cpuid);
   %mus(ff,:) = mu(ii);
    %for i=1:n
    %    plot(time(ff),mus(i),'.');
    %end

   %getframe(figM);
   %grid on;
   fclose(fid);
   
   ff = ff+1;
    
end
 





%hold on;

fig(1);clf;
    
set(gcf,'Color',[1 1 1]);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,40,20]);
set(gcf,'PaperPositionMode','auto');

for i=1:n
    plot(shear*time(:),mus(:,i));
    axis([0 max(shear*time) 0 2000]);
    set(gca,'TickLength',[0.02 0.02],'FontSize',40);
    ylabel('$\mu$'); xlabel('$St$');
    set(gca,'Units','normalized','Position',[0.15,0.15,0.82,0.82]);
    plotTickLatex2D('ytickdy',0.005,'ytickdx',0.003,'xtickdy',-0.015,'ylabeldx',-0.02,'xlabeldy',-0.02);
    %pause(1);

    saveas(gcf,['/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/track/png/img',numlab],'png');
end



%{
[nn,xout] = hist(log10(0.5*vsq),100);
xlog = 10.^xout;
ylog = nn.*sqrt(xlog);
loglog(xlog,ylog/max(max(ylog)));
axis([0.1 1e4 1e-4 2]);
%}
%plot(x1,v1,'.');

%{
xx = 0:8:2048;
vv = -100:10:100;
mvx = zeros(n,2);
for i=1:n
    mvx(i,1) = x1(i);
    mvx(i,2) = v1(i);
end
%}

%nvx = hist2d(mvx,vv,xx);

%
%   nXBins = length(vXEdge);
%   nYBins = length(vYEdge);
%   vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
%   vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
%   pcolor(vXLabel, vYLabel,mHist2d); colorbar



%plot(x1,x2,'.');
%axis([x1l x1u x2l x2u]);
%axis equal;
% view(0,0);
% grid on;

%{
fig(1);clf;
[nn,xout] = hist(vprl,1000);
plot(xout,nn);
%bar(xout,nn,'barwidth',1,'basevalue',1)
set(gca,'YScale','log')

fig(2);clf;
[nn,xout] = hist(v2,1000);
plot(xout,nn);
%bar(xout,nn,'barwidth',1,'basevalue',1)

fig(3);clf;
[nn,xout] = hist(v3,1000);
plot(xout,nn);
%bar(xout,nn,'barwidth',1,'basevalue',1)

fig(4);clf;
xrange = min(v2):0.1:max(v2);
yrange = min(v1):0.1:max(v1);
data(:,1) = v2(:);
data(:,2) = v1(:);

    for i=1:length(xrange)-1

        data((data(:,1)>xrange(i))&(data(:,1)<=xrange(i+1)),3)=i;
    end

    for i=1:length(yrange)-1

        data((data(:,2)>yrange(i))&(data(:,2)<=yrange(i+1)),4)=i;  

    end

    count=zeros(length(xrange)-1,length(yrange)-1);

    data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
    data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away
    
    for i=1:size(data,1)
        count(data(i,3),data(i,4))=count(data(i,3),data(i,4))+1; 
    end

imagesc(xrange,yrange,count);
axis xy; axis image; shading flat;

%}