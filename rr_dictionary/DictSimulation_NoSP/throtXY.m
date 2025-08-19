function Rth=throtXY(phi,theta,dir)

Rz = zrot(-theta);
if dir == 0
    Rx = xrot(phi);
    Rth = inv(Rz)*Rx*Rz;
else
    Ry = yrot(phi);
    Rth = inv(Rz)*Ry*Rz;
end
