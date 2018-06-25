function calc_scattering_vpar()
basename='./';
basename='/tigress/dstonge/pegasus/sat_b4/track/';
files=dir(sprintf('%s*.bin',basename));
npart = length(files(:));

mt=0;
MAX_T=500;

framesize = 10;
iframe = 1.0/framesize;

smoothing=1;

NBINS=ceil(MAX_T/framesize);
data = zeros(1,NBINS);
data2 = zeros(1,NBINS);
counts = zeros(1,NBINS);
firsts = zeros(1,npart);

for i = 1:npart
  name=files(i).name;
  fp = fopen(sprintf('%s%s',basename,name),'r');
  x = fread(fp,1,'int');
  y = fread(fp,1,'int');
  A = fread(fp,[x y],'single');
  fclose(fp);
  
  time = A(1,:); 
  vpar = A(9,:); 

  ind=0; 
  ind0=0; 
  curdt=0;
  t0 = 0;
  j=1;
  vp0=vpar(1);
  cont=1;

  nsteps = length(time(:));
  while ( j < nsteps)
    ind = floor(time(j)*iframe);
    if(ind > mt) mt = ind; end;
    if(ind ~= ind0)
      if(cont > 0 && cont < 0)
        curdt = time(j) - t0;
        for k =(j+1):nsteps
          curdt = curdt + (time(k) - time(k-1));
          if(vpar(k)*vp0 < 0)
            break;
          end
        end
        data(ind) = data(ind) + curdt;
        data2(ind) = data2(ind) + 1.0/curdt;
        if(ind0 == 0) 
          firsts(i) = curdt;
        end
      end
      ind0 = ind;
      t0 = time(j);
      vp0 = vpar(j);
      cont=1;
    else
      if(cont > 0 && (vpar(j)*vp0 < 0))
        curdt = time(j) - t0;
        if(ind0==0)
          firsts(i) = curdt;
        end
        data(ind+1)= data(ind+1) + curdt;
        data2(ind+1)= data2(ind+1) + 1.0/curdt;
        counts(ind+1)=counts(ind+1) + 1;
        cont = -1;
      end
    end
    j=j+1;
  end
  if(mod(i,100) == 0) 
    disp(i);
  end
end


fp = fopen(sprintf('%sscatt_vpar',basename),'w');
for i = 1:mt
  fprintf(fp,'%e %d %e %e\n',(i-0.5)*framesize,counts(i),npart/data(i), data2(i)/npart);
end
fclose(fp);

fp = fopen(sprintf('%svpar_firsts',basename),'w');
for i=1:npart
  fprintf(fp,'%d %e\n',i,firsts(i));
end
fclose(fp);

