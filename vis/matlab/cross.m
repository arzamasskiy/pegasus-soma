function [cx,cy,cz] = cross(ax,ay,az,bx,by,bz)
    cx = ay.*bz - az.*by;
    cy = az.*bx - ax.*bz;
    cz = ax.*by - ay.*bx;
end
