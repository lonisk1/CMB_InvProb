function [A] = getRayTraceMat2(x_s, y_s, x_d, y_d, rad, T, E, n, n_E, n_g, visual, scale)
%
% [A] = getRayTraceMat2(x_s, y_s, x_d, y_d, rad, T, E, n, n_E, n_g, visual, scale)
%
%
% This function generates a sparse matrix A that represents cone-beam
% ray-tracing.  Bilinear interpolation is used to determine the weights in
% the matrix.  Here we have a higher resolution true image (2x the spatial
% directions)
%
%   Output: A - sparse matrix
%   Input:  x_s, y_s - source dimensions
%           x_d, y_d - detector dimensions
%           rad, T   - radius of circle and height of detector
%           E        - height of event planes
%           n, n_E   - 3D grid size
%           n_g      - initial grid
%           visual   - 0 to only save the observed detectors
%                      1 to save the entire grid of detectors (for
%                      visualization purposes)
%           scale    - scale for the image grid (e.g., 2)
%

A = []; % initialization
A_horzLayer = []; % initialization

% Scale the detector and source grid x_d,y_d
x_s = scale*x_s; y_s = scale*y_s;
x_d = scale*x_d; y_d = scale*y_d;
xvec_i = -max(x_s(:)) : max(x_s(:));
yvec_i = -max(y_s(:)) : max(y_s(:));
[x_i,y_i] = ndgrid(xvec_i, yvec_i);
n_g = scale*n_g;
rad = scale*rad;

for j = 1:size(x_s,1) % loop over one dimension of source
  for k = 1:size(y_s,2) % loop over second dimension of source
    xs = x_s(j,1);
    ys = y_s(1,k);

    % Identify all detectors using circle with radius
    idx = find((x_d-xs).^2 + (y_d-ys).^2 <= (rad+(scale)*0.5)^2 & (x_d-xs).^2 + (y_d-ys).^2 >= (rad-(scale)*0.5)^2);

    %% For each of the detectors identified in idx, determine location of intersection at slice located at E_j
    ii_row = []; jj_row = []; ss_row = [];
    for u = 1:n_E % loop over event planes
      ii = []; jj = []; ss = []; % re-initialize each time to form 'local' spare matrix
      for i = 1:length(idx) % loop over detector pixels
        % for each detector, find the point of intersection at event time E
        intx = xs + (E(u)/T * x_d(idx(i))); % x-value of ray through event plane from source xs
        inty = ys + (E(u)/T * y_d(idx(i))); % y-value of ray through event plane from source ys

        if abs(intx) <= n_g && abs(inty) <= n_g % check if there is an intersection with the ray
          intx_lo = floor(intx); % lower x-integer
          inty_lo = floor(inty); % lower y-integer
          intx_up = ceil(intx); % upper x-integer
          inty_up = ceil(inty); % upper y-integer

          % 4 different cases for ray intersection with event plane
          if intx == intx_lo && inty == inty_lo
            %warning('ray passed through a single node...all is lost')
            % in this case all the weight should be passed to the
            % corresponding node
            a1 = 1;
            j1 = find(x_i == intx_lo & y_i == inty_lo);
            if visual
              ii = [ii; idx(i)];
            else
              ii = [ii; i];
            end
            jj = [jj; j1];
            ss = [ss; a1];

          elseif intx == intx_lo %(or intx_up: it doesn't matter which)
            % in this case there should only be two weights as the ray passes
            % between two nodes who have the same x-coord. (use line weighting)
            % Diagram:
            %                 (x_int,y_int)
            %                     |
            %                     v a1
            %                     -----
            %     (x_lo,y_lo) X---*---X (x_lo,y_up)
            %                 -----
            %                   a2
            a1 = inty_up - inty;
            a2 = 1 - a1;
            j1 = find(x_i == intx_lo & y_i == inty_lo);
            j2 = find(x_i == intx_lo & y_i == inty_up);

            if visual
              ii = [ii; idx(i); idx(i)];
            else
              ii = [ii; i; i];
            end
            jj = [jj; j1; j2];
            ss = [ss; a1; a2];

          elseif inty == inty_lo %(or inty_up: it doesn't matter which)
            % in this case there should only be two weights as the ray passes
            % between two nodes who have the same y-coord. (use line weighting)
            %   Diagram:
            %        {  X (x_lo,y_lo)
            %   a2   {  |
            %        {  * (x_int,y_int)   }
            %           |                 } a1
            %           X (x_up,y_lo)     }
            %
            a1 = intx_up - intx;
            a2 = 1 - a1;
            j1 = find(x_i == intx_lo & y_i == inty_lo);
            j2 = find(x_i == intx_up & y_i == inty_lo);

            if visual
              ii = [ii; idx(i); idx(i)];
            else
              ii = [ii; i; i];
            end
            jj = [jj; j1; j2];
            ss = [ss; a1; a2];

          else
            % in this case there should be four weights as the ray passes
            % between four nodes (use bilinear interpolation)
            %
            % Diagram:
            % (x_lo,y_lo)X---------X (x_lo,y_up)    * <-- (x_int,y_int)
            %            | a2 | a3 |
            %            |----*----|
            %            | a1 | a4 |
            % (x_up,y_lo)X---------X (x_up,y_up)
            %
            % Bilinear interp. - finding weights "a_i"
            a1 = (intx_up - intx) * (inty - inty_lo); % corresponds to j2
            a2 = (intx - intx_lo) * (inty - inty_lo); % corresponds to j4
            a3 = (intx - intx_lo) * (inty_up - inty); % corresponds to j3
            a4 = (intx_up - intx) * (inty_up - inty); % corresponds to j1


            j1 = find(x_i == intx_lo & y_i == inty_lo);
            j2 = find(x_i == intx_lo & y_i == inty_up);
            j3 = find(x_i == intx_up & y_i == inty_lo);
            j4 = find(x_i == intx_up & y_i == inty_up);
            if visual
              ii = [ii; idx(i); idx(i); idx(i); idx(i)];
            else
              ii = [ii; i; i; i; i];
            end
            jj = [jj; j1; j2; j3; j4];
            ss = [ss; a4; a1; a3; a2]; %obscure order is because of sketch in notes
          end % 4 cases
        end % intersection with ray
      end % loop over idx

      ii_row = [ii_row; ii];
      jj_row = [jj_row; (u-1)*length(x_i(:)) + jj ];
      ss_row = [ss_row; ss];

    end % loop over event planes
    if ~isempty(ii_row)
      if visual
        A_horzLayer = sparse(ii_row, jj_row, ss_row, size(x_d,1)*size(x_d,2), size(x_i,1)*size(x_i,2)*n_E);
      else
        A_horzLayer = sparse(ii_row, jj_row, ss_row, length(idx), size(x_i,1)*size(x_i,2)*n_E);
      end
      A = cat(1, A, A_horzLayer); %concatinate vertically as a function of source location
      A_horzLayer = []; %re-initialize
    end
  end % end source second dimension
end % end source first dimension