function spec_to_growth()
basename='/home/dstonge/PegData/prod_run3';
files=dir(sprintf('%s/B_*.dat',basename));
nspec = length(files(:));
dt=4;
t=0;
A = dlmread(sprintf('%s/%s',basename,files(1).name),'')';
ll = length(A(2,:));
k=A(1,:);
norm = A(2,2);
data=zeros(nspec,ll+1);


for i = 1:nspec
   name=files(i).name;
   A = dlmread(sprintf('%s/%s',basename,name),'')';
   
   data(i,1) = t;
   data(i,2:(ll+1)) = A(2,:)/norm;
   t = t+dt;
end



fp=fopen(sprintf('%s/modal_B.dat',basename),'w');
for i = 1:nspec
    for j = 1:(ll+1)
        fprintf(fp,'%e  ', data(i,j));
    end
    fprintf(fp,'\n');
end

fclose(fp);

ldata = log(data(:,2:(ll+1)));
ddata = ((circshift(ldata,-1) - ldata  )/dt)';
ddata = ddata(:,(1:nspec-1));
dtdata = zeros(length(k),nspec);
dtdata(:,1) = k;
dtdata(:,2:nspec) = ddata;

fp=fopen(sprintf('%s/modal_B_gr.dat',basename),'w');
for i = 1:ll
    for j = 1:nspec
        fprintf(fp,'%e  ', dtdata(i,j));
    end
    fprintf(fp,'\n');
end

fclose(fp);



files=dir(sprintf('%s/V_*.dat',basename));
nspec = length(files(:));
t=0;
A = dlmread(sprintf('%s/%s',basename,files(1).name),'')';
ll = length(A(2,:));
norm = A(2,2);
k=A(1,:);
data=zeros(nspec,ll+1);

for i = 1:nspec
   name=files(i).name;
   A = dlmread(sprintf('%s/%s',basename,name),'')';
   
   data(i,1) = t;
   data(i,2:(ll+1)) = A(2,:)/norm;
   t = t+dt;
end


fp=fopen(sprintf('%s/modal_V.dat',basename),'w');
for i = 1:nspec
    for j = 1:(ll+1)
        fprintf(fp,'%e  ', data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);


ldata = log(data(:,2:(ll+1)));
ddata = ((circshift(ldata,-1) - ldata  )/dt)';
ddata = ddata(:,(1:nspec-1));
dtdata = zeros(length(k),nspec);
dtdata(:,1) = k;
dtdata(:,2:nspec) = ddata;


fp=fopen(sprintf('%s/modal_V_gr.dat',basename),'w');
for i = 1:ll
    for j = 1:nspec
        fprintf(fp,'%e  ', dtdata(i,j));
    end
    fprintf(fp,'\n');
end

fclose(fp);