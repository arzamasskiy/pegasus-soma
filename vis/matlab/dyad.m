function [c11,c12,c13,c21,c22,c23,c31,c32,c33] = dyad(ax,ay,az,bx,by,bz)
    c11= ax.*bx;
    c12= ax.*by; 
    c13= ax.*bz;
    c21= ay.*bx;
    c22= ay.*by; 
    c23= ay.*bz;
    c31= az.*bx;
    c32= az.*by; 
    c33= az.*bz;
end
