function ascii_to_bin()
basename='./';
basename='/tigress/dstonge/pegasus/sat_b4_proto_3/track/';
files=dir(sprintf('%sdata/*.lis',basename))


npart = length(files(:))
name=files(1).name;

for i = 1:npart
  name=files(i).name;
  A = dlmread(sprintf('%sdata/%s',basename,name),'\t',3,0)';

  d=size(A);

  fp = fopen(sprintf('%sdata/turb.%06d.bin',basename,i-1),'w');
  fwrite(fp,d(1),'int');
  fwrite(fp,d(2),'int');
  fwrite(fp,single(A),'single');
  fclose(fp);
  
  if(mod(i,100) == 0) 
      disp(i);
  end
end
