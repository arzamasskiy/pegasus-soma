function calc_scattering()
basename='./';
files=dir(sprintf('%s*.lis',basename));

MAX_T=500;

%binsize = 0.1;
%ibin = 1.0/binsize;

framesize = 2;
iframe = 1.0/framesize;

mindt=0.01;
logsize=1.25;

smoothing=11;


thres=3.0;
diff_thres=0.075;

NBINS=45;
TBINS=floor(MAX_T/framesize)+2;

data = zeros(NBINS,TBINS);

centers = zeros(2,TBINS-1);
raw = zeros(2,TBINS-1);
raw_d = zeros(2,TBINS-1);

centers(1,:) = 0:framesize:MAX_T;
raw(1,:) = 0:framesize:MAX_T;
raw_d(1,:) = 0:framesize:MAX_T;
npart = length(files(:));
%data(:,1) = binsize*(1:NBINS);
data(:,1) = 0.01*logsize.^(0:(NBINS-1));


for i = 1:100
    name=files(i).name;
    A = dlmread(sprintf('%s%s',basename,name),'\t',3,0);
  
    time= A(:,1);
    mu = smooth(A(:,8),smoothing); 
   
  
    deri = smooth(diff(mu),smoothing);
   
    nsteps = length(time(:));
    dlmwrite(sprintf('%sderi_%d',basename,i),[time(2:nsteps) mu(2:nsteps) deri],' ');
    
    deri_scal = deri./mu(2:nsteps);
    
    mu0=mu(1);
    d0=0;
    t0 = time(1);
    for j = 1:nsteps
        t1i = floor(iframe*time(j)) + 2;
        if( mu0 > (thres*mu(j)) || mu(j) > thres*mu0)
            dt = time(j)-t0;
            
            t0i = floor(iframe*t0) + 2; % 1 to go from starting index 0 to 1, another to add axis info
            %bi = min(floor(dt*ibin)+1,NBINS); %linear
            bi = min(max(floor(log(dt/mindt)/log(logsize)),0)+1,NBINS); %logarithmic

            raw(2,t1i) = raw(2,t1i) + 1;
            
            for k = t0i:t1i
              data(bi,k) = data(bi,k) + 1;
            end
            t0 = time(j);
            mu0 = mu(j);     
        end
        if(j > 1)
            if(abs(deri_scal(j-1)) > diff_thres && d0 <= diff_thres)
                raw_d(2,t1i) = raw_d(2,t1i) + 1;
            end
            d0 = deri_scal(j-1);
        end
    end
    if(mod(i,100) == 0) 
        disp(i);
    end
end

for i = 2:TBINS
 f=fit((1:NBINS).',data(:,i),'gauss1');
 centers(2,i-1) = mindt*logsize^(f.b1-1);
end

raw(2,:) = raw(2,:)./(npart*framesize);
raw_d(2,:) = raw_d(2,:)./(npart*framesize);

dlmwrite(sprintf('%shistogram',basename),data ,' ');
dlmwrite(sprintf('%stimes',basename),centers' ,' ');
dlmwrite(sprintf('%srawtimes',basename),raw' ,' ');
dlmwrite(sprintf('%srawtimes_d',basename),raw_d' ,' ');