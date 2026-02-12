<<<<<<< HEAD
% rotate(x,y,z) by agl (in radius)
function [nx, ny, nz] = CoordRotate(x,y,z,agl)
[tx,ty,p2xy] = CoordTrans(x,y,z);
ny = ty;
nx = tx*cos(p2xy+agl);
nz = tx*sin(p2xy+agl);
end

% transfer Cartesian coordinate system to rotation plane system
function [tx,ty,p2xy] = CoordTrans(x,y,z)
    ty = y;
    if x~=0
    tx = sqrt(x^2+z^2)*(x/abs(x));
    p2xy = atan(z/x); % angle plane to x-y palne in radius; + --> clockwise
    else
        if z~=0
            tx = z;
            p2xy = pi/2;
        else
            tx = 0;
            p2xy = 0;
        end
    end
=======
% rotate(x,y,z) by agl (in radius)
function [nx, ny, nz] = CoordRotate(x,y,z,agl)
[tx,ty,p2xy] = CoordTrans(x,y,z);
ny = ty;
nx = tx*cos(p2xy+agl);
nz = tx*sin(p2xy+agl);
end

% transfer Cartesian coordinate system to rotation plane system
function [tx,ty,p2xy] = CoordTrans(x,y,z)
    ty = y;
    if x~=0
    tx = sqrt(x^2+z^2)*(x/abs(x));
    p2xy = atan(z/x); % angle plane to x-y palne in radius; + --> clockwise
    else
        if z~=0
            tx = z;
            p2xy = pi/2;
        else
            tx = 0;
            p2xy = 0;
        end
    end
>>>>>>> 5e4cf7b (Initial commit)
end