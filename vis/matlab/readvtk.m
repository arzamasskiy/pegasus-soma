function [time,var,name,status] = readvtk(Grid,filename,varid)
% 
% readvtk:  READS A .vtk FILE AND RETURNS THE DESIRED VARIABLE.

% OPEN FILE FOR READING
[fid, message] = fopen(filename,'r');
if (fid==-1)
    fprintf(2,'[readvtk]:  %s could not be opened!\n', filename);
    fprintf(2,'%s', message);
    status = -1;
    return;
end;

% READ HEADER
fgets(fid);
fscanf(fid,'%s',4);
time = fscanf(fid,'%e',1);

fgets(fid);fgets(fid);fgets(fid);
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

fscanf(fid,'%s',1);
fscanf(fid,'%d',1);

% CHECK FOR CONSISTENCY
warn = false;
warn = warn && (nx1   == Grid.nx1);
warn = warn && (nx2   == Grid.nx2);
warn = warn && (nx3   == Grid.nx3);
warn = warn && (dx1   == Grid.dx1);
warn = warn && (dx2   == Grid.dx2);
warn = warn && (dx3   == Grid.dx3);
warn = warn && (time  == Grid.time  );
warn = warn && (x1min == Grid.x1min );
warn = warn && (x2min == Grid.x2min );
warn = warn && (x3min == Grid.x3min );
if (warn)
  fprintf(2,'[readvtk]:  %s failed consistency check!\n',filename);
  status = -1;
  return;
end;

% SET THE FILE POINTER TO THE BEGINNING OF THE DATA
for i=1:varid-1
 if (Grid.ary(i) == 1) % scalar
   typ = fscanf(fid,'%s',1);
   name = fscanf(fid,'%s',1);
   if (~strcmp(typ,'SCALARS'))
    fprintf(2,'[readvtk]:  %s is not a scalar!\n', name);
    status = -1;
    return;
   end
   fgets(fid);fgets(fid);
   offset = Grid.ndata*4;
   fseek(fid,offset,'cof');
 end
 if (Grid.ary(i) == 3) % vector
   typ = fscanf(fid,'%s',1);
   name = fscanf(fid,'%s',1);
   if (~strcmp(typ,'VECTORS'))
    fprintf(2,'[readvtk]:  %s is not a vector!\n', name);
    status = -1;
    return;
   end
   fgets(fid);
   offset = 3*Grid.ndata*4;
   fseek(fid,offset,'cof');
 end
 if (Grid.ary(i) == 9) % tensor
   typ = fscanf(fid,'%s',1);
   name = fscanf(fid,'%s',1);
   if (~strcmp(typ,'TENSORS'))
    fprintf(2,'[readvtk]:  %s is not a tensor!\n', name);
    status = -1;
    return;
   end
   fgets(fid);
   offset = 9*Grid.ndata*4;
   fseek(fid,offset,'cof');
 end
end

% READ DATA
i=varid;
if (Grid.ary(i) == 1) % scalar
  typ = fscanf(fid,'%s',1);
  name = fscanf(fid,'%s',1);
  if (~strcmp(typ,'SCALARS'))
   fprintf(2,'[readvtk]:  %s is not a scalar!\n', name);
   status = -1;
   return;
  else
    fgets(fid);fgets(fid);
    var = reshape(fread(fid,Grid.ndata,'float','b'),nx1,nx2,nx3);
    fprintf('[readvtk]:  %s read, first value = %f\n',name,var(1,1));
  end
end
if (Grid.ary(i) == 3) % vector
  typ = fscanf(fid,'%s',1);
  name = fscanf(fid,'%s',1);
  if (~strcmp(typ,'VECTORS'))
    fprintf(2,'[readvtk]:  %s is not a vector!\n', name);
    status = -1;
    return;
  else
    fgets(fid);
    var = reshape(fread(fid,3*Grid.ndata,'float','b'),3,nx1,nx2,nx3);
    fprintf('[readvtk]:  %s read, first value = %f\n',name,var(1,1));
  end
end
if (Grid.ary(i) == 9) % tensor
  typ = fscanf(fid,'%s',1);
  name = fscanf(fid,'%s',1);
  if (~strcmp(typ,'TENSORS'))
    fprintf(2,'[readvtk]:  %s is not a tensor!\n', name);
    status = -1;
    return;
  else
    fgets(fid);
    var = reshape(fread(fid,9*Grid.ndata,'float','b'),9,nx1,nx2,nx3);
    fprintf('[readvtk]:  %s read, first value = %f\n',name,var(1,1));
  end
end

% CLOSE FILE
status = fclose(fid);
if (status == -1)
    fprintf(2,'[readvtk]:  %s could not be closed!\n', filename);
end;

return;