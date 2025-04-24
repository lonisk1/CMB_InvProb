function xtrue = getstring(n, E, params)
%
% xtrue = getstring(m,n,E,type, params)
%
% This function constructs a 3D image that is n x n x E based on the
% function in params.fn
%   
% params is a structure giving additional parameters
%
% J. Chung and L. Onisk, 4/2025

example = params.example;
c = params.c;
a = params.a;

% get grids
xgrid = linspace(-3,3,n);
ygrid = linspace(-3,3,n);
tgrid = linspace(0,1,E);

% get full space
[x_obj, y_obj, z_obj] = ndgrid(xgrid, ygrid, tgrid); % 3D object grid

switch example
    case 1
        xtrue = double(1 - sqrt((x_obj - c*z_obj).^2 + (y_obj).^2));
        xtrue (xtrue < 0) = 0;
        xtrue = xtrue.^a;
    case 2
        xtrue = double(1 + .1* sin(2*pi*z_obj) - sqrt((x_obj).^2 + (y_obj).^2)); % change 0.1 to 0.5
        xtrue (xtrue < 0) = 0;
        xtrue = xtrue.^a;
    case 3
        z_obj = 4*z_obj; % scale z to [0,4]
        xtrue1 = double(1 - .5*(z_obj) - sqrt((x_obj).^2 + (y_obj).^2));
        xtrue1 (xtrue1 < 0) = 0;
        xtrue2 = double(1 - .5*(4-z_obj) - sqrt((x_obj).^2 + (y_obj).^2));
        xtrue2 (xtrue2 < 0) = 0;
        xtrue = zeros(n,n,E);
        xtrue(:,:,1:floor(E/2)) = xtrue1(:,:,1:floor(E/2));
        xtrue(:,:,floor(E/2)+1:end) = xtrue2(:,:,floor(E/2)+1:end);
        xtrue = xtrue.^a;

  case 4 % binary sphere example
        tgrid = linspace(0,1,E-2*lessSlices);
        % get adapted space
        [x_obj, y_obj, z_obj] = ndgrid(xgrid, ygrid, tgrid); % 3D object grid

        radius = 0.5;
        xc = 0; yc = 0; zc = 0.5; % the center of sphere
        logicalSphere = (x_obj-xc).^2 + (y_obj-yc).^2 + (z_obj-zc).^2 <=radius*radius;
        xtrue = zeros(n,n,E-2*lessSlices);
        xtrue(logicalSphere) = 1; % set to zero

        filler = zeros(n,n,lessSlices);
        xtrue = cat(3,filler,xtrue); %concatinate front-end
        xtrue = cat(3,xtrue,filler); %concatinate back-end

  case 5 % smoothed sphere example
        xgrid = linspace(-3,3,n); 
        ygrid = linspace(-3,3,n);
        tgrid = linspace(0,2.5,E-2*lessSlices);
         % get adapted space
        [x_obj, y_obj, z_obj] = ndgrid(xgrid, ygrid, tgrid); % 3D object grid

        xc = 0; yc = 0; zc = 1.25; % the center of sphere
        xtrue = double(1 - sqrt((x_obj).^2 + (y_obj).^2 + (z_obj-zc).^2));
        xtrue (xtrue < 0) = 0;
        xtrue = xtrue.^a;

        filler = zeros(n,n,lessSlices);
        xtrue = cat(3,filler,xtrue); %concatinate front-end
        xtrue = cat(3,xtrue,filler); %concatinate back-end
        
  case 6 % floating sphere
        xgrid = linspace(-3,3,n); 
        ygrid = linspace(-3,3,n);
        tgrid = linspace(0,4,E);
         
        % get adapted space
       [x_obj, y_obj, z_obj] = ndgrid(xgrid, ygrid, tgrid); % 3D object grid

        xc = 0; yc = 0; zc = 2; % the center of sphere
        rad = 1.5;
        xtrue = double(rad - sqrt((x_obj-xc).^2 + (y_obj-xc).^2 + (z_obj-zc).^2));
        xtrue (xtrue < 0) = 0;
        xtrue = xtrue.^a;
    otherwise
        error('Need to set example in params.')
end

