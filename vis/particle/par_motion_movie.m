clear all;
figlab = 3;
figM = figure(figlab); hold off;
fname='mrshear'; % file name
fn1=2000;    % first figure index
fn2=6280;  % last figure index

Lx = 1152;
Ly = 1152;
shear = 3e-4;
dir ='/Volumes/My Passport Studio/pegasus/mrs-slow/vtk/';        % directory
fname = 'mrshear';

tracked = 1;
it=337;
% 1032, 1814, 2272, 127, 6, 74, 257, 258, 268, 271, 278, 289, 304, 327, 328,
% 330, 394, 408, 556, 620, 861, 899, 916, 1111, 1126, 1621, 1733

%high = 0;
%low  = 0;

for fp=fn1:fn2
    
    fp
 if (fp<10)
     numlab = ['000',num2str(fp)];
     else if (fp<100)
         numlab = ['00',num2str(fp)];
         else if (fp<1000)
             numlab = ['0',num2str(fp)];
             else
                 numlab = num2str(fp);
             end
         end
 end

 
 %fid=fopen(['/Users/kunz/Documents/codes/pegasus/bin/',fname,'.',numlab,'.par.lis'],'rb');
 %fid =fopen(['/Users/kunz/track/',fname,'.',numlab,'.track.lis'],'rb');
 fid =fopen(['/Volumes/My Passport Studio/pegasus/mrs-slow/track/',fname,'.',numlab,'.track.lis'],'rb');
 
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
     x1(j) = parinfo(1); x2(j) = parinfo(2); %x3(j) = parinfo(3);
     %v1(j) = parinfo(4); v2(j) = parinfo(5); v3(j) = parinfo(6);
     %f_0(j) = parinfo(7); mu(j) = parinfo(8);
     %vprl(j) = parinfo(9); vprp(j) = parinfo(10);
     prop(j) = fread(fid,1,'int32');
     pid(j) = fread(fid,1,'int64');
     cpuid(j) = fread(fid,1,'int64');

 end
 
 [cpu,ii] = sort(cpuid);
 
 x1p = x1(ii);
 x2p = x2(ii);
 
 %{
 
 for j=1:npar
     dB(ff-fn1*10+1,j) = interp2(Grid.x1,Grid.x2,dBprl',x1(ii(j)),x2(ii(j)));
 end
 %}
 
 
 %mus(ff,:) = mu(ii);
 
 %clear x1 x2 x3 v1 v2 v3 f_0 mu vprl vprp prop pid cpuid;

 
 plot(x1p(it),x2p(it),'.'); hold off;
 axis([x1l x1u x2l x2u]);
 getframe(figM);
 %view(0,0);
 grid on;
 
 fclose(fid);
 
 %getframe(figM);
 
 %view(0,0);
% grid on;
% hold on;
 
 end

