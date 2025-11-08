function Z = zernike_noll(j, rho, theta)
% 前8项Zernike（未归一化，Noll编号）
switch j
    case 1,  Z = ones(size(rho));                          % piston
    case 2,  Z = rho.*cos(theta);                          % tilt x
    case 3,  Z = rho.*sin(theta);                          % tilt y
    case 4,  Z = 2*rho.^2 - 1;                             % defocus
    case 5,  Z = rho.^2.*sin(2*theta);                     % astig 45°
    case 6,  Z = rho.^2.*cos(2*theta);                     % astig 0/90°
    case 7,  Z = (3*rho.^3-2*rho).*cos(theta);             % coma x
    case 8,  Z = (3*rho.^3-2*rho).*sin(theta);             % coma y
    otherwise, error('Zernike项数量哆啦');
end
end
