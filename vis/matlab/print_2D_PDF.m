function [hist,xedges,yedges]=print_2D_PDF(file,datx,daty,logx,logy)
  if (length(datx(:))~=length(daty(:)))
    fprintf(2, 'List size mismatch \n');
    return;
  end

  inx = datx;
  iny=daty;
  if(logx)
     inx=log10(datx); 
  end
  if(logy)
     iny=log10(daty); 
  end
     
  [hist,xedges,yedges]=histcounts2(inx,iny);
  lx=length(hist(:,1));
  ly=length(hist(1,:));
  xwidth=xedges(2)-xedges(1);
  ywidth=yedges(2)-yedges(1);
  
  mass=sum(hist(:));
  bEx= xedges(1:lx) + xwidth*0.5; 
  bEy= yedges(1:ly) + ywidth*0.5; 
  
  vals = hist/mass;
  
  
  if(logx)
     xwidths=10.^(xedges(2:(lx+1)))-10.^(xedges(1:lx));
     vals = vals./repmat(xwidths',1,ly);
     %bEx=10.^bEx; gnuplot doesn't like these
  else
     vals = vals/xwidth; 
  end
  if(logy)
     ywidths=10.^(yedges(2:(ly+1)))-10.^(yedges(1:ly));
     vals = vals./repmat(ywidths,lx,1);
     %bEy=10.^bEy; gnuplot doesn't like these
  else
    vals = vals/ywidth;
  end
  
  maxHist = 1.0/max(vals(:));
  NVal=vals*maxHist;
  
  %NtVal = NVal.*(NVal >= thres);
  
  histout = fopen(file,'W');
  for j = 1:ly
    for i = 1:lx
       fprintf(histout,'%e %e %e %e\n',bEx(i),bEy(j), ...
                                        vals(i,j), NVal(i,j));
    end
    fprintf(histout,'\n');
  end
  fclose(histout);
end
