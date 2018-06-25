function calc_scattering_smooth()
basename='./';
basename='/tigress/dstonge/pegasus/sat_b4_proto_3/track/';
files=dir(sprintf('%sdata/*.bin',basename));

om = dlmread(sprintf('%som.dat',basename),' ',0,0)';
tsmooth=pi./om(2,:);

% 

MAX_T=700;
thres=exp(1);
thres=1.25;

framesize = 20;
iframe = 1.0/framesize;

NBINS=MAX_T/framesize;
data = zeros(1,NBINS);
data2 = zeros(1,NBINS);
counts = zeros(1,NBINS);
firsts = zeros(1,NBINS);


name=files(1).name;
fp = fopen(sprintf('%sdata/%s',basename,name),'r');
x = fread(fp,1,'int');
y = fread(fp,1,'int');
A = fread(fp,[x y],'single')';
fclose(fp);
time= A(:,1);
nsteps = length(time(:));

mu_s = zeros(1,nsteps);

clear b;
n=1;
mt=0;


npart = length(files(:));
npart = 2000;

vals = zeros(NBINS,npart);

for i = 1:npart
    name=files(i).name;
    fp = fopen(sprintf('%sdata/%s',basename,name),'r');
    x = fread(fp,1,'int');
    y = fread(fp,1,'int');
    A = fread(fp,[x y],'single')';
    fclose(fp);

%   mu = smooth(A(:,8),smoothing,'lowess'); 
    mu = (A(:,8)); 
    time= A(:,1);
   

    %smooth over gyroperiod?
    s_l=1;
    s_r=1;
    p=1;

    for j=1:nsteps
      t=time(j);
      dt=tsmooth(j);
      lt=t-dt;
      rt=t+dt;
      while(time(s_l) < lt)
        s_l=s_l+1;
      end
      while(s_r < nsteps & time(s_r+1) < rt)
        s_r=s_r+1;
      end
      s_s=0;
      s_s1=0;
      ind= max(1,(s_r - j) - (j-s_l));
      mult=1;
      for k=s_l:s_r
        s_s = s_s+ ind*log10(mu(k));
        s_s1=s_s1+ind;
        if(k ==j) mult=-1; end;
        ind=ind+mult;
      end
      %fprintf('%d %d %d\n',j-s_l,s_r-j,s_s1);
      %mu_s(j) = 10^(s_s/s_s1);
       mu_s(j) = exp((sum(log(mu(s_l:s_r)).^p)/(s_r-s_l+1))^(1/p));

    end
    
    if(i<5)
      fp = fopen(sprintf('%stester_%d',basename,i),'w');
      for j = 1:nsteps
        fprintf(fp,'%e %e %e %e\n',time(j),mu_s(j),mu(j),mu(j)-mu_s(j));
      end
      fclose(fp);
    end
    mu=mu_s;

    ind=0; 
    ind0=0; 
    curdt=0;
    t0 = 0;
    j=1;
    mu0=mu(1);

    while ( j < nsteps)
      ind = floor(time(j)*iframe);
      if(ind > mt) mt = ind; end;
      if(ind ~= ind0)
        if(mu0 > 0)
          curdt = time(j) - t0;
          for k =(j+1):nsteps
            curdt = curdt + (time(k) - time(k-1));
            if k == nsteps
              sprintf('oof! %s\n',files(i).name)
            end
            if(mu(k)/mu0 > thres || mu(k)/mu0 < 1/thres)
              break;
            end
          end
          data(ind) = data(ind) + curdt;
          data2(ind) = data2(ind) + 1.0/curdt;
          vals(ind,i) = curdt;
        end
        ind0 = ind;
        t0 = time(j);
        mu0 = mu(j);
      else
        if(mu0 > 0 && (mu(j)/mu0 > thres || mu(j)/mu0 < 1/thres))
          curdt = time(j) - t0;
          data(ind+1)= data(ind+1) + curdt;
          data2(ind+1)= data2(ind+1) + 1.0/curdt;
          vals(ind+1,i) = curdt;
          counts(ind+1)=counts(ind+1) + 1;
          mu0 = -1;
        end
      end
      j=j+1;
    end

    if(mod(i,100) == 0) 
        disp(i);
    end
end

fp = fopen(sprintf('%sscatt_mu_om_test',basename),'w');
fprintf(fp,'0 0 0\n');
for i = 1:mt
  fprintf(fp,'%e %d %e %e\n',(i-0.5)*framesize,counts(i),npart/data(i), data2(i)/npart);
end
fclose(fp);


for i= 1:mt
  [nv1,edges1] = histcounts(vals(i,:), 'Normalization', 'pdf');
  [nv2,edges2] = histcounts(1.0./vals(i,:),'Normalization','pdf');
  l1=length(nv1(:));
  l2=length(nv2(:));

  fp = fopen(sprintf('%shist_t%d',basename,framesize*i),'w');
  for j = 1:l1
    fprintf(fp,'%e %e\n',0.5*(edges1(j)+edges1(j+1)),nv1(j));
  end
  fclose(fp);

  fp = fopen(sprintf('%shistinv_t%d',basename,framesize*i),'w');
  for j = 1:l2
    fprintf(fp,'%e %e\n',0.5*(edges2(j)+edges2(j+1)),nv2(j));
  end
  fclose(fp);

end
