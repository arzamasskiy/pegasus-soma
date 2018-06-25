function [Grid,status] = init_grid(filename,ary)
% 
% init_grid:  INITIALIZES A GRID STRUCTURE BY OPENING A .vtk FILE AND
% READING IN HEADER INFORMATION ONLY. SOME EXTRA VARIABLES ARE COMPUTED
% FOR CLARITY.

% OPEN FILE FOR READING
[fid, message] = fopen(filename,'r');
if (fid==-1)
    fprintf(2,'[init_grid]:  %s could not be opened!\n', filename);
    fprintf(2,'%s', message);
    status = -1;
    return;
end;

% READ HEADER
fgets(fid);
fscanf(fid,'%s',4);
time = fscanf(fid,'%e',1);

fgets(fid);
fgets(fid);
fgets(fid);
fscanf(fid,'%s',1);
nx1 = fscanf(fid,'%d',1)-1;
nx2 = fscanf(fid,'%d',1)-1;
nx3 = fscanf(fid,'%d',1)-1;
ndim = (nx1>1)+(nx2>1)+(nx3>1);
if ndim <= 2
    nx3 = 1;
end
if ndim == 1
    nx2 = 1;
end

fscanf(fid,'%s',1);
x1min = fscanf(fid,'%e',1);
x2min = fscanf(fid,'%e',1);
x3min = fscanf(fid,'%e',1);

fscanf(fid,'%s',1);
dx1 = fscanf(fid,'%e',1);
dx2 = fscanf(fid,'%e',1);
dx3 = fscanf(fid,'%e',1);

x1max = x1min + nx1*dx1;
x2max = x2min + nx2*dx2;
x3max = x3min + nx3*dx3;

fscanf(fid,'%s',1);
ndata = fscanf(fid,'%d',1);

if (ndata ~= nx1*nx2*nx3)
  fprintf(2,'[init_grid]:  data size inconsistent!\n');
  status = -1;
  return;
end;

% CLOSE FILE
status = fclose(fid);
if (status == -1)
    fprintf(2,'[init_grid]:  %s could not be closed!\n', filename);
end;

% VARIABLE ARRANGEMENT
nvar=max(size(ary));
for i=1:nvar
  if ((ary(i) ~= 1) && (ary(i) ~= 3) && (ary(i) ~= 9))
    fprintf(2, '[init_grid]:  variable arrangement inconsistent!\n');
    status = -1;
    return;
  end
end

% COMPUTE SOME DERIVED QUANTITIES
%x1i = linspace(x1min,x1max,nx1+1);
%x2i = linspace(x2min,x2max,nx2+1);
%x3i = linspace(x3min,x3max,nx3+1);
x1 = linspace(x1min+0.5*dx1,x1max-0.5*dx1,nx1);
x2 = linspace(x2min+0.5*dx2,x2max-0.5*dx2,nx2);
x3 = linspace(x3min+0.5*dx3,x3max-0.5*dx3,nx3);

% INITIALIZE THE REST OF GRID STRUCTURE
Grid.nx1         = nx1;
Grid.nx2         = nx2;
Grid.nx3         = nx3;
Grid.dx1         = dx1;
Grid.dx2         = dx2;
Grid.dx3         = dx3;
Grid.x1min       = x1min;
Grid.x1max       = x1max;
Grid.x2min       = x2min;
Grid.x2max       = x2max;
Grid.x3min       = x3min;
Grid.x3max       = x3max;
Grid.ndim        = ndim;
Grid.nvar        = nvar;
Grid.ndata       = ndata;
Grid.ary         = ary;
Grid.x1          = x1;
Grid.x2          = x2;
Grid.x3          = x3;
Grid.time        = time;

status = 0;

return;