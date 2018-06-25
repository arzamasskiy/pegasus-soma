nbrcores = 512;  % specify number of processors
prlbins = linspace(-4.0,4.0,400);
% construct parallel-velocity bins (see header in *.spec files)
prpbins = linspace(0.0,4.0,200);
% construct perpendicular-velocity bins
nbr_cells = 400*200;  % total number of bins

% march through processors
for numlab=1:nbrcores-1
 nstr=num2str(numlab);
file=['/tigress/dstonge/pegasus/run7_alt2_spec/id',nstr,'/dynamo-id',nstr,'.spec'];
% file name
fid=fopen(file,'r');
% open file
C = textscan(fid,repmat('%f',[1,1+nbr_cells]),'CommentStyle','#');
% read file, ignoring header
fclose(fid);
% close file
D = horzcat(C{:});
% reform data into 2d (vprl,vprp) matrix
[m,n] = size(D);
f = D(m,2:n); 
f = reshape(f,[400,200]); 
imagesc(prlbins,prpbins,f');
axis image; axis xy;
title(['Processor =',nstr]);
drawnow;
end
