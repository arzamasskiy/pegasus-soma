function calc_Bgrowthrates()
  basename='./';
  A = dlmread(sprintf('%sdynamo_clean.hst',basename),'',3,0)';
  time = A(1,:);
  bgr = smooth(log(A(11,:) + A(12,:) + A(13,:)),11);
  
  
  %plot(time,bgr);
  %drawnow;
   
  l = length(time(:));
  d_bgr = zeros(l,1);
  
  for i=2:l
     d_bgr(i) = 0.5*(bgr(i)-bgr(i-1))/(time(i)- time(i-1)); 
  end
  plot(time,d_bgr);
  drawnow
  
  dlmwrite(sprintf('%sgrowth_rates.dat',basename),[time' d_bgr] ,' ');
end