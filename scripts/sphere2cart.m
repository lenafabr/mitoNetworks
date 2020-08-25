function cartpts = sphere2cart(spherepts,cent)
% convert from spherical coordinates (r, theta, phi) to cartesian
% sphere at center cent

cartpts = spherepts(:,1).*[sin(spherepts(:,2)).*cos(spherepts(:,3)), ...
    sin(spherepts(:,2)).*sin(spherepts(:,3)),cos(spherepts(:,2))] + cent;
end