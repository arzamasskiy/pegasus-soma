function calc_scattering_mu()
basename='./';
basename='/tigress/dstonge/pegasus/prod_run3/track/';
files=dir(sprintf('%sdata/*.bin',basename))

MAX_T=480;
thres=exp(1);

framesize = 16;
iframe = 1.0/framesize;

smoothing=1;


NBINS=MAX_T/framesize;
data = zeros(1,NBINS);
data2 = zeros(1,NBINS);
counts = zeros(1,NBINS);
firsts = zeros(1,NBINS);



clear b;
n=1;
mt=0;

npart = length(files(:));

vals = zeros(NBINS,npart);

for i = 1:npart
    name=files(i).name;
    fp = fopen(sprintf('%sdata/%s',basename,name),'r');
    x = fread(fp,1,'int');
    y = fread(fp,1,'int');
    A = fread(fp,[x y],'single')';
    fclose(fp);

    %mu = smooth(A(:,8),smoothing,'lowess'); 
    mu = A(:,8); 
    time= A(:,1);
   
    ind=0; 
    ind0=0; 
    curdt=0;
    t0 = 0;
    j=1;
    mu0=mu(1);

    nsteps = length(time(:));
    while ( j < nsteps)
      ind = floor(time(j)*iframe);
      if(ind > mt) mt = ind; end;
      if(ind ~= ind0)
        if(mu0 > 0)
          curdt = time(j) - t0;
          for k =(j+1):nsteps
            curdt = curdt + (time(k) - time(k-1));
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

fp = fopen(sprintf('%sscatt_mu_tcorr',basename),'w');
fprintf(fp,'0 0 0\n');
for i = 1:mt
  fprintf(fp,'%e %d %e %e\n',(i-0.5)*framesize,counts(i),npart/data(i), data2(i)/npart);
end
fclose(fp);

for i= 1:mt
  [nv1,edges1] = histcounts(vals(i,:), 'Normalization', 'pdf');
  l1=length(nv1(:));

  fp = fopen(sprintf('%shist_t%d',basename,framesize*i),'w');
  for j = 1:l1
    fprintf(fp,'%e %e\n',0.5*(edges1(j)+edges1(j+1)),nv1(j));
  end
  fclose(fp);

end
