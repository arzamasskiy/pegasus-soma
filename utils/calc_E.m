function calc_Bgrowthrates()
  basename='./';
  basename='/Users/dstonge/PegData/s_small_run2_redo/';
  A = dlmread(sprintf('%sdynamo_clean.hst',basename),'',3,0)';
  time = A(1,:);
  %bgr = smooth(log(A(11,:) + A(12,:) + A(13,:)),11);
  TE = 0.5*(2.0*A(21,:) + A(20,:));
  KE = (A(8,:) + A(9,:) + A(10,:));
  ME = (A(11,:) + A(12,:) + A(13,:));
  E = TE+KE+ME;
  %plot(time,bgr);
  %drawnow;
   
  l = length(time(:));
  dTE = zeros(l,1);
  dME = zeros(l,1);
  dKE = zeros(l,1);
  dE = zeros(l,1);
  
  for i=2:l
     dTE(i)= (TE(i)-TE(i-1))/(time(i)- time(i-1)); 
     dME(i)= (ME(i)-ME(i-1))/(time(i)- time(i-1)); 
     dKE(i)= (KE(i)-KE(i-1))/(time(i)- time(i-1)); 
     dE(i) = (E(i)-E(i-1))/(time(i)- time(i-1)); 
  end
  
  fp = fopen(sprintf('%sdEdt.dat',basename),'w');
  fprintf(fp,'[1] t [2] E [3] dEdt [4] dTEdt [5] dKEdt [6] dMEdt\n');   
  for i=1:l
    fprintf(fp,'%e %e %e %e %e %e\n',time(i),E(i),dE(i),dTE(i),dKE(i),dME(i));   
  end
  fclose(fp);
end