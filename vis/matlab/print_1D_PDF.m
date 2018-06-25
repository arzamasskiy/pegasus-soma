function [hist,edges]=print_1D_PDF(file,dat,logx)
  in=dat;
  if(logx)
     in=log10(dat);
  end
    
  [hist,edges]=histcounts(in);
  l=length(hist(:));
  width=edges(2)-edges(1);
  
  mass=sum(hist(:));
  bE= edges(1:l) + width*0.5; 
  
  vals = hist/mass;
  
  if(logx)
     vals = vals./(10.^(edges(2:(l+1)))-10.^(edges(1:l)));
     bE=10.^bE;
  else
     vals = vals/width;
  end
  
  maxHist = 1.0/(max(vals(:)));
  NVal=vals*maxHist;
  
  histout = fopen(file,'W');

  for ii = 1:l
    fprintf(histout,'%e %e %e\n',bE(ii),vals(ii), NVal(ii));
  end
  fclose(histout);
end

