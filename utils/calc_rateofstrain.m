function calc_rateofstrain()
basename='./';
files=dir(sprintf('%sV_spec_*.dat',basename));


nfiles = length(files(:))';

time = 4*(0:(nfiles-1))';
ros =  zeros(1,nfiles);

A = dlmread(sprintf('%sV_spec_0000.dat',basename),' ')';
noise = 2*A(2,:).*A(1,:).^2;

l=length(noise);
n = zeros(1,l);
for i = 1:l
    n(i) = n(max(i-1,1))+noise(i);
end
plot(n)
drawnow

for i = 1:nfiles
    name=files(i).name;
    A = dlmread(sprintf('%s%s',basename,name),' ')';
    l=length(A);
    B = zeros(3,l);
       
    
    for j = 1:l
       B(1,j) = A(1,j); 
       B(2,j) = B(2,max(j-1,1)) + A(2,j);
       B(3,j) = B(3,max(j-1,1)) + 2*A(2,j)*A(1,j)^2;
    end
    B(3,:) = 0.36*sqrt((B(3,:)-n));
    ros(i) = B(3,l);
    dlmwrite(sprintf('%s%s.new',basename,name),B',' ');
end


plot(time,ros)
drawnow;
  dlmwrite(sprintf('%sROS.dat',basename),[time ros'] ,' ');
  