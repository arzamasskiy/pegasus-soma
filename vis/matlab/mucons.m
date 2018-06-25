clear all;
figlab = 4;
figM = figure(figlab); hold off;
fn1=1;    % first figure index
fn2=6280;  % last figure index

Lx = 1152;
Ly = 1152;
shear = 3e-4;
dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/track/';%'/Volumes/My Passport Studio/pegasus/mrs-slow/track/';        % directory
fname = 'mrshear';

tracked = -1;

%high = 0;
%low  = 0;

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
 
 fid=fopen(['/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/track/',fname,'.',numlab,'.track.lis'],'rb');
 %fid =fopen(['/Users/kunz/track/',fname,'.',numlab,'.track.lis'],'rb');
 %fid =fopen([dir,fname,'.',numlab,'.track.lis'],'rb');
 
 % Read the coordinate limits
 coorlim = fread(fid,12,'float');
 x1l = coorlim(1); x1u = coorlim(2); x2l = coorlim(3); x2u = coorlim(4); %x3l = coorlim(5); x3u = coorlim(6);
 %x1dl = coorlim(7); x1du = coorlim(8); x2dl = coorlim(9); x2du = coorlim(10); x3dl = coorlim(11); x3du = coorlim(12);
 
 % Read number of particle types
 npartypes = fread(fid,1,'int32');
 for i=1:npartypes
     m(i) = fread(fid,1,'float');
     qomc(i) = fread(fid,1,'float');
 end
 
 % Read the time
 time = fread(fid,1,'float');
 dt   = fread(fid,1,'float');
 
 % Read the particle number
 npar =  fread(fid,1,'int64');
 
 % Read all the particle information
 for j=1:npar
     parinfo = fread(fid,10,'float');
     x1(j) = parinfo(1); x2(j) = parinfo(2);% x3(j) = parinfo(3);
     v1(j) = parinfo(4); v2(j) = parinfo(5); v3(j) = parinfo(6);
     %f_0(j) = parinfo(7); 
     mu(j) = parinfo(8);
     vprl(j) = parinfo(9); 
     vprp(j) = parinfo(10);
     prop = fread(fid,1,'int32');
     pid = fread(fid,1,'int64');
     cpuid(j) = fread(fid,1,'int64');
     %if (cpuid(j)==tracked)  
     %    if (x1(j) >= Grid.x1(kx) && x1(j) < Grid.x1(kx+1))
     %       if (x2(j) >= Grid.x2(ky) && x2(j) < Grid.x2(ky+1))
         
         %xc = floor(x1(j))+0.5*Lx+1;
         %yc = floor(x2(j))+0.5*Ly+1;
         %dB = dBprl(xc,yc)

     %    dB(ff-fn1*10+1) = interp2(Grid.x1,Grid.x2,dBprl',x1(j),x2(j));

     %        end
     %    end
         %}
         %plot(x1(j),x2(j),'.');
         %axis([x1l x1u x2l x2u]);
         %xc(ff-fn1+1) = x1(j);
         %yc(ff-fn1+1) = x2(j);
     %end
 end
 
 
 [cpu,ii] = sort(cpuid);
 %{
 for j=1:npar
     dB(ff-fn1*10+1,j) = interp2(Grid.x1,Grid.x2,dBprl',x1(ii(j)),x2(ii(j)));
 end
 %}
 
 
 x1s(ff,:) = x1(ii);
 x2s(ff,:) = x2(ii);
 v1s(ff,:) = v1(ii);
 v2s(ff,:) = v2(ii);
 v3s(ff,:) = v3(ii);
 mus(ff,:) = mu(ii);
 vprls(ff,:) = vprl(ii);
 vprps(ff,:) = vprp(ii);
 
 
 
 
 %plot(x1,x2,'.');
 %axis equal;
 %hold on;
 %axis([x1l x1u x2l x2u]);
 
 
 clear x1 x2 v1 v2 v3 mu vprl vprp cpuid;

 %getframe(figM);
 
 %view(0,0);
% grid on;
% hold on;
 
 fclose(fid);
 
end


B = (vprps.^2) ./ mus;



    shear = 3e-4;
    dir ='/Users/kunz/Documents/codes/pegasus/bin/mrs-slow/';        % directory
    fname = 'mrshear';
    filename = [dir,fname,'.hst'];

    fid = fopen(filename);
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
    fclose(fid);

    time = C{1};
    Bxa  = C{17};
    Bya  = C{18};
    
    
    B0 = sqrt(Bxa.^2+Bya.^2);
    Bm = B0(2:6281);
    
    for i=1:2304
        dB(:,i) = B(:,i) - Bm(:);
    end
    
    
    
    
    
    