% Thomas C. Day
% Generating a snowflake yeast, or many, according to these parameters.
% This script is particularly meant for Kai's ploidy project.
% Vary aspect ratio, and vary cell size separately.

%% Inputs:
% folder = 'C:\Users\Thomas\Dropbox (GaTech)\Yunker Lab\Data\Project Ploidy\';
% cd(folder);

V               = 225; %4/3*pi*(2.5^3); % this is the volume of every cell
N               = 1; % number of clusters to generate at each parameter value
numGens         = 20; % number of generations of cell division in a group, here i pick a high value
AR              = 1%:.1:2; % aspect ratio value
err_AR          = 0; % standard deviation in aspect ratio
diam            = 2 .* ((3*V./(4*pi*AR)).^(1/3)); % smallest diameter of the cell, scaled by AR
err_diam        = 0; % variation in cell size
pole_theta      = deg2rad(10); % buds nearest the pole will be chosen from between 0 and 10 degrees in polar angle
THETA           = deg2rad(45); % polar angle average from SEM data
thetaVariance   = deg2rad(0); % variation in polar angle from SEM data
distance_thresh = 1.1672 * (diam/4.58); % minimum distance (um) separating bud scars, this grows proportionately to cell diameter
new_bud_prob    = .8; % probability that the first cell will bud near the pole
check_overlap   = 1; % do we check the overlap?
overlap_thresh  = 1e2 * ((diam(1)/4.58).^2); % threshold of total overlaps, this scales with cell volume but not with aspect ratio

% Do we visualize figures?
figure_viz = 0;

% Do we save to file?
para_save = 0;

% Initialization:
fprintf('Welcome to the cell simulator v2.0\n');
fprintf('Written by Thomas C. Day, 2020\n');

% Save parameters as a mat file:
% foldername = ['CellVol=',num2str(V,3),'\'];
% mkdir(foldername);
% save([foldername,'parameters.mat']);

%% MAIN CODE: -------------------------------------------------------------
% -------------------------------------------------------------------------
% Parallel iteration over varied params:
for a = 1:length(AR) % sweeping cell aspect ratio
    cells_sim = cell(N,1);
    aspRat    = AR(a);
    theta     = THETA;
    filename = ['thomas-sim','_AspRat=',num2str(aspRat,'%1.1f'),...
                '_CellDiam=',num2str(diam(a),3),...
                '_pole-theta=',num2str(rad2deg(pole_theta),3),'_theta=',num2str(rad2deg(theta),3),'_dtheta=',num2str(rad2deg(thetaVariance),3),...
                '_back-prob=',num2str(new_bud_prob,3),...
                '_N=',num2str(N)];
    filename = strrep(filename, '.','-');
    COM = zeros(3,N);
    RAD = zeros(1,N);

    % Loop over N clusters
    for ii = 1:N
        
        % Generate one cluster:
        fprintf([num2str(ii), ' / ', num2str(N), '\n']);
        [cell_list] = ELYES_SIM(diam(a), err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh(a), new_bud_prob, numGens, 0, check_overlap, overlap_thresh);
               
        % Add newest cell list to main list
        cells_sim{ii} = cell_list; % add to list
        
        % Save data to file:
        if para_save == 1
            if mod(ii,N) == 0 % save every N clusters generated
                parsave([foldername,filename,'.mat'], cells_sim);
            end
        end
        
    end
end

% Show figure of one cluster:
if figure_viz == 1
    figure;
    hold on; box on; set(gca,'linewidth',2);
    for k = 1:length(cell_list)
        [x,y,z] = VISUALIZE_ELLIPSOID(cell_list(k), 30);
        surf(x,y,z,'facealpha',1,'edgecolor','none');
    end
    view(3); axis equal;
    lighting gouraud;
    lightangle(0,30);
    material dull;
end

%% FUNCTIONS
function [x,y,z] = VISUALIZE_ELLIPSOID(cell_of_interest, resolution)
%%C: the covariance matrix.
%%Dir: direction of the estimate, to be plotted together with the DT ellipsoid.
%%M: the mean vector, usually 0 in case of DT.
%%speed: time to pause between plotting, lowervalue = faster.
%%Dir: 1 or -1 for clockwise, anticlockwise.
%%time: number of iterations for which rotation is needed, higher = longer.
%%example: visualizeDTrot(diag([17 2 2]),[0 0 0],0.4,1,100)

% Extract info from cell_of_interest:
radii = cell_of_interest.Radii;
centers = cell_of_interest.Center;
R = cell_of_interest.Rmatrix;

% Generate data for "unrotated" ellipsoid
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3), resolution);

% Rotate data with orientation matrix R and center T
a = kron(R(:,1), xc);
b = kron(R(:,2), yc);
c = kron(R(:,3), zc);
data = a+b+c; n = size(data,2);

% Store for output:
x = data(1:n,:) + centers(1); 
y = data(n+1:2*n,:) + centers(2); 
z = data(2*n+1:end,:) + centers(3);

end

function [cell_list] = ELYES_SIM(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, numGens, neighbor_thresh, check_overlap, overlap_thresh)

    cell_list       = [];  % this will eventually be the tabulated list
    
    % Create root cell:
    rootCell.Center       = [0; 0; 0];
    rootCell.Radii        = 1/2 * [aspRat * diam; diam; diam];
    rootCell.Rmatrix      = eye(3);
    rootCell.Generation   = 0;
    rootCell.IDnumber     = 0;
    rootCell.Daughters    = [];
    rootCell.Mother       = [];
    rootCell.BudXYZ       = [];
    rootCell.Overlaps     = [];
    rootCell.Neighborhood = [];
    cell_list             = [cell_list, rootCell];
    counting_ix           = 0; % this will count how many cells are in the cluster

    % Iterate over generations:
    for g = 1:numGens
        
        % First check the total amount of overlap in the cluster:
        fprintf(['Gen = ',num2str(g),'\n']);
        if check_overlap == 1
            flag = CHECK_OVERLAPS2(cell_list, overlap_thresh);
        else
            flag = 0;
        end
        
        if flag ~= 1
            % If there is not too much overlap:
            nPossible   = length(cell_list); % number of possible new cells
            variedTheta = theta + thetaVariance * randn(nPossible, 1); % list of theta values corresponding to these cells

            % Assign positions for the newest cells:
            for n = 1:nPossible
                [newBud, newBudRel, newAxis] = getDaughterPos(cell_list(n), pole_theta, variedTheta(n), distance_thresh, new_bud_prob); % finds the new bud xyz position and relative position to the old cell
                if newBud ~= 0 % if the budding chance was successful
                    counting_ix = counting_ix + 1; % the number of cells goes up by one

                    % record new bud in mother cell information array:
                    cell_list(n).BudXYZ = [cell_list(n).BudXYZ, newBudRel];
                    cell_list(n).Daughters = [cell_list(n).Daughters, counting_ix];

                    % record new cell information for the cell list:
                    newCell.IDnumber     = counting_ix;
                    newCell.Generation   = g;
                    newCell.Mother       = cell_list(n).IDnumber;
                    newCell.Daughters    = [];
                    newCell.BudXYZ       = [];
                    newCell.Overlaps     = [];
                    newCell.Neighborhood = [];
                    newCell.Radii        = 1/2 * [aspRat * diam; diam; diam];

                    % Find new cell center:
                    [Sm, Rm, Tm]   = GET_SURFACE_MATRICES(cell_list(n), .1);
                    A              = newCell.Radii(1);
                    p              = A * newAxis + newBudRel;
                    q              = Tm * Rm * [p; 1];
                    newCell.Center = q(1:3);

                    % Choose a cell orientation based upon surface normal axis:
                    avec   = newCell.Center - newBud(1:3);
                    avec   = avec./norm(avec);
                    bvec_x = rand;
                    bvec_y = rand;
                    bvec_z = - (avec(1)*bvec_x + avec(2)*bvec_y)/avec(3);
                    bvec   = [bvec_x; bvec_y; bvec_z];
                    bvec   = bvec./norm(bvec);
                    cvec   = cross(avec, bvec);
                    newCell.Rmatrix = [avec, bvec, cvec];

                    % add new cell to cell list:
                    cell_list = [cell_list, newCell];
                end
            end
        else
            cell_list = cell_list;
            break;
        end
        

    end
    
end

function [flag] = CHECK_OVERLAPS2(cell_list, overlap_thresh)

    % Obtain distances between all particles:
    [D, ~, ~, ~] = get_particle_distances(cell_list);
    
    % Obtain combined equatorial radii of two particles:
    radii = [cell_list.Radii];
    eq_radii = radii(2,:);
    R = eq_radii + eq_radii';
    
    % Obtain the energy associated with each interaction:
    d = D - R;
    d = d .* ~eye(size(d));
    w = d(d < 0);
    u = w.^2;
    U = sum(u,'all'); % sum up all interaction energies
    
    flag = U > overlap_thresh; % if U is too big, flag is thrown
    
end

function [D, X, Y, Z] = get_particle_distances(cell_list)
    % Determine distances between each pair of particles:
    r_in = [cell_list.Center]';
    X = r_in(:,1) - r_in(:,1)';
    Y = r_in(:,2) - r_in(:,2)';
    Z = r_in(:,3) - r_in(:,3)';
    D = sqrt(X.^2 + Y.^2 + Z.^2);
end

function [flag] = CHECK_OVERLAPS(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude, overlap_thresh)

    [~, ~, ~, ~, overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude);
    Overlaps = sum(overlaps, 2);
    if ~isempty(Overlaps)
        flag = Overlaps(2) > overlap_thresh;
    else
        flag = 0;
    end
    

end

function parsave(filename, cell_list)
    save(filename, 'cell_list');
end

function [S, R, T] = GET_SURFACE_MATRICES(cell_of_interest, scaling)

% Scaling matrix
S = [cell_of_interest.Radii(1) + scaling, 0, 0, 0;...
        0, cell_of_interest.Radii(2) + scaling, 0, 0; ...
        0, 0, cell_of_interest.Radii(3) + scaling, 0; ...
		0,0,0,1];
%  cell wall is usually ~75 nm thick, this accounts for that

% rotation matrix
R = cell_of_interest.Rmatrix;
R = [R; 0,0,0]; R = horzcat(R,[0;0;0;1]);

% translation matrix
T = eye(4); 
T(1,end) = cell_of_interest.Center(1);
T(2,end) = cell_of_interest.Center(2);
T(3,end) = cell_of_interest.Center(3);

end

function [newPos, relativePos, axis] = getDaughterPos(CELL, POLE_TH, TH, DISTANCE_THRESH, NEW_BUD_PROB_THRESH)

% Surface definition for possible mother cell:
[Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);

% Determine where on cell body to bud the next scar:
if isempty(CELL.Daughters) % reproduce near the pole with prob 0.7
    frac = rand;
    if frac < NEW_BUD_PROB_THRESH
        % bud at pole
        th = POLE_TH * rand;
        ph = 2 * pi * rand;
        ix_too_close = [];
        r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
        x = r * cos(th);
        y = r * sin(th) * cos(ph);
        z = r * sin(th) * sin(ph);
    else
        % bud on the polar angle
        ph = 2 * pi * rand;
        th = TH;
        ix_too_close = [];
        r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
        x = r * cos(th);
        y = r * sin(th) * cos(ph);
        z = r * sin(th) * sin(ph);
    end

else % reproduce based upon the polar angle theta, sometimes back bud (if 3 or more bud scars)
    
    nDaughters = size(CELL.BudXYZ,2);
    if nDaughters < 4
        ph = 2 * pi * rand;
        th = TH;
    else
        frac = rand;
        if frac < 1 % bud at the distal pole
            ph = 2 * pi * rand;
            th = TH;
        else % bud at proximal pole
            ph = 2 * pi * rand;
            th = deg2rad(180) - TH;
        end
    end

    % Check if new location is too close to existing scars:
    r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
    x = r * cos(th);
    y = r * sin(th) * cos(ph);
    z = r * sin(th) * sin(ph);
    t = [x; y; z];
    s = CELL.BudXYZ - t;
    d = sqrt(s(1,:).^2 + s(2,:).^2 + s(3,:).^2);
    ix_too_close = find(d < DISTANCE_THRESH);
end

if isempty(ix_too_close) % new bud is successful
    relativePos = [x; y; z];
    newPos      = Tm * Rm * [relativePos; 1];
    normal      = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
    axis        = normal./norm(normal);
else % this bud is too close, loses its chance
    relativePos = 0;
    newPos      = 0;
    axis        = 0;
end

end

function [cell_list_in, cell_list_out, Overlaps] = INCLUDE_FORCES(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude, mobility_pos, mobility_rot, dt, T, figure_viz)
% This function takes as input a grown snowflake and calculates the forces
% and torques on each cell. From these forces, it the calculates any
% rearrangements.

cell_list_in = cell_list;

for t = 1:T
    fprintf(['Time = ',num2str(t),'\n']);
    [Forces, Torques, flog1, flog2, Overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude);
    [cell_list] = UPDATE_POSITIONS(cell_list, Forces, Torques, mobility_pos, mobility_rot, dt);
    
    if figure_viz == 1
        if mod(t,1) == 0
            close all;
            figure('visible','on','units','centimeters','position',[1,1,18.3,18.3]); 
            hold on; box on; set(gca,'linewidth',2);
            for n = 1:length(cell_list)
                c_o_i = cell_list(n);
                [x,y,z] = VISUALIZE_ELLIPSOID(c_o_i, 30);
                surf(x,y,z,'facealpha',1,'edgecolor','none');
                %{
%                 Fnet = 3*Forces(:,n) + c_o_i.Center;
%                 Tnet = 3*Torques(:,n) + c_o_i.Center;
%                 plot3([c_o_i.Center(1), Fnet(1)], [c_o_i.Center(2),Fnet(2)], [c_o_i.Center(3), Fnet(3)],'r-','linewidth',3);
%                 plot3([c_o_i.Center(1), Tnet(1)], [c_o_i.Center(2),Tnet(2)], [c_o_i.Center(3), Tnet(3)],'b-','linewidth',3);
                forces = flog1{n};
                ix = find(vecnorm(forces) ~= 0);
                forces(:,ix) = forces(:,ix)./vecnorm(forces(:,ix));
                pts = flog2{n};
                plot3(pts(1,1), pts(2,1), pts(3,1), 'rx','markersize',12,'linewidth',2); % plot birth scar location
                plot3([pts(1,1), pts(1,1)+forces(1,1)], [pts(2,1), pts(2,1)+forces(2,1)], [pts(3,1), pts(3,1)+forces(3,1)],'r-','linewidth',2); % plot force from birth scar
                nDaughters = length(cell_list(n).Daughters);
                if nDaughters == 0
                    nDaughters = 1;
                end
                for jj = 1:nDaughters
                    plot3(pts(1,jj+1), pts(2,jj+1), pts(3,jj+1),'b.','markersize',10,'linewidth',2);
                    plot3([pts(1,jj+1), pts(1,jj+1)+forces(1,jj+1)], [pts(2,jj+1), pts(2,jj+1)+forces(2,jj+1)], [pts(3,jj+1), pts(3,jj+1)+forces(3,jj+1)],'b-','linewidth',2); % plot force from daughters
                end
                plot3(pts(1,nDaughters+2:end), pts(2,nDaughters+2:end), pts(3,nDaughters+2:end),'go','markersize',14,'linewidth',2);
                plot3([pts(1,nDaughters+2:end), pts(1,nDaughters+2:end) + forces(1,nDaughters+2:end)],[pts(2,nDaughters+2:end), pts(2,nDaughters+2:end) + forces(2,nDaughters+2:end)], [pts(3,nDaughters+2:end), pts(3,nDaughters+2:end) + forces(3,nDaughters+2:end)],'g-','linewidth',2);
                %}
            end
            view(3); axis equal;
            lightangle(15,15);
            lighting gouraud;
            material dull;
            centers = [cell_list.Center];
            xlim([min(centers(1,:)) - 5, max(centers(1,:)) + 5]);
            ylim([min(centers(2,:)) - 5, max(centers(2,:)) + 5]);
            zlim([min(centers(3,:)) - 5, max(centers(3,:)) + 5]);
            print(['test_sims_t=',num2str(t,'%03.f')],'-dpng','-r500');
        end
    end
end
    
cell_list_out = cell_list;

end


function [cell_list_out] = UPDATE_POSITIONS(cell_list, Forces, Torques, mobility_pos, mobility_rot, dt)
% Overdamped motion
cell_list_out = cell_list;
dx = mobility_pos * Forces * dt; % magnitude and direction of displacement
zero_magnitude = zeros(size(vecnorm(Torques)));
da = - mobility_rot * vecnorm(Torques) * dt; % magnitude of rotational change
Tvec = Torques./vecnorm(Torques); % axis of rotation

% Update cell positions:
for n = 1:length(cell_list)
    cell_list_out(n).Center = cell_list(n).Center + dx(:,n);
end

% Update cell orientations:
for n = 1:length(cell_list)
    if da(n) ~= 0 % only rotate if there is a torque
        Rrot = rotation_around_arb_axis(Tvec(:,n), da(n));
        cell_list_out(n).Rmatrix = Rrot * cell_list(n).Rmatrix;
    end
end

end


function [Forces, Torques, forces_log1, forces_log2, final_overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude)

% Obtain distances between all particles:
[D, ~, ~, ~] = get_particle_distances(cell_list);

% Find all pairs of cells within the neighborhood threshold:
N = D < neighbor_thresh; % finds all neighbor cells
N = N - eye(size(N)); % we don't care about counting a cell as its own neighbor

% Obtain forces/torques for each cell:
Forces = zeros(3, length(cell_list));
Torques = zeros(3, length(cell_list));
final_overlaps = [];
for n = 1:length(cell_list)
    neighbors = find(N(:,n) == 1); % only need to check the neighboring cells for interactions
    
    % Steric interactions:
    [F_pts, O_pts, Overlaps, Directions] = get_overlaps(cell_list, n, neighbors);
    final_overlaps = [final_overlaps, Overlaps];
    
    % Chitin interactions:
    [m_pts, d_pts] = CHITIN_INTERACTIONS(cell_list, n);
    
    % Calculate forces:
    % Force from chitin bond with mother cell
    if isempty(m_pts)
        F_mother = [0;0;0];
    else
        F_mother = chitin_magnitude * m_pts;
    end
    % Forces from chitin bonds with daughters
    if isempty(d_pts)
        F_daughter = [0;0;0];
    else
        F_daughter = chitin_magnitude * d_pts;
    end
    % Forces from steric interactions
    if ~isempty(Overlaps)
        F_steric = steric_magnitude * Overlaps(2,:) .* Directions;
    else
        F_steric = [0;0;0];
    end
    % All forces
    forces = [F_mother, F_daughter, F_steric];
    Forces(:,n) = sum(forces, 2);
    
    % Log forces and force pts for visualization:
    forces_log1{n} = forces;
    [S,R,T] = GET_SURFACE_MATRICES(cell_list(n),0);
    M = T*R*S;
    mother_log = M*[-1;0;0;1];
    if ~isempty(cell_list(n).Daughters)
        daughter_log = T*R*[cell_list(n).BudXYZ; ones(1,length(cell_list(n).Daughters))];
        if isempty(F_pts)
            F_pts = mother_log(1:3);            
        end
        forces_log2{n} = [mother_log(1:3), daughter_log(1:3,:), F_pts];
    else
        if isempty(F_pts)
            F_pts = mother_log(1:3);
        end
        forces_log2{n} = [mother_log(1:3), mother_log(1:3), F_pts];
    end
    
    % Calculate torques:
    % First calculate the r-vector:
    [S, R, T]  = GET_SURFACE_MATRICES(cell_list(n), 0); 
    M = T*R*S;
    if isempty(cell_list(n).Mother)
        r_mother = [0;0;0;1];
    else
        r_mother   = M * [-1; 0; 0; 1] - M * [0;0;0;1];
    end
    if ~isempty(cell_list(n).Daughters)
        r_daughter = T * R * [cell_list(n).BudXYZ; ones(1,size(cell_list(n).BudXYZ,2))] - M*[0;0;0;1];
    else
        r_daughter = [0;0;0;1];
    end
    r_steric   = [F_pts; ones(1,size(F_pts,2))] - M*[0;0;0;1];
    if isempty(r_steric)
        r_steric = [0;0;0;1];
    end
    r_angle = -r_mother(1:3);
    if ~isempty(cell_list(n).Mother)
        cell_mother = cell_list(n).Mother + 1;
        ix = find(cell_list(cell_mother).Daughters == n - 1);
        budxyz = cell_list(cell_mother).BudXYZ(:,ix);
        [s, r, t] = GET_SURFACE_MATRICES(cell_list(cell_mother),0);
        surf_norm = 2 * [budxyz(1)/s(1,1)^2; budxyz(2)/s(2,2)^2; budxyz(3)/s(3,3)^2];
        preferred_axis = r * [surf_norm; 1];
        current_axis = R * S * [1;0;0;1];
        preferred_axis = preferred_axis(1:3)/norm(preferred_axis(1:3));
        current_axis = current_axis(1:3)/norm(current_axis(1:3));
        v = bond_torque_magnitude * (preferred_axis - current_axis);
    else
        v = [0;0;0];
    end
    forces = [forces, v];
    rdisp = [r_mother(1:3,:), r_daughter(1:3,:), r_steric(1:3,:), r_angle];
    torques = zeros(size(forces));
    for i = 1:size(forces,2)
        torques(:,i) = cross(forces(:,i), rdisp(:,i));
    end
    Torques(:,n) = sum(torques, 2);
end

end

function [Force_pts, Overlap_pts, Overlaps, Directions] = get_overlaps(cell_list, ix_o_i, nbors)
    
    % Obtain overlapping volume and center point of overlap:
    cell1       = cell_list(ix_o_i);
    Force_pts   = [];
    Overlap_pts = [];
    Overlaps    = [];
    Directions  = [];
    for j = 1:length(nbors)
        cell2 = cell_list(nbors(j));
        
        % Check if the two cells overlap at all: --------------------------  
        % First numerically approximate the surface of each cell
        [x1,y1,z1] = VISUALIZE_ELLIPSOID(cell1, 20); % numerical surface
        [x2,y2,z2] = VISUALIZE_ELLIPSOID(cell2, 20);
        x1 = x1(:); y1 = y1(:); z1 = z1(:);
        x2 = x2(:); y2 = y2(:); z2 = z2(:);
        r1 = [x1,y1,z1,ones(size(x1))]';
        r2 = [x2,y2,z2,ones(size(x2))]';
        
        % Second find the analytic surface of each ellipsoid
        [S1,R1,T1] = GET_SURFACE_TRANSFORMATIONS_SIM(cell1,.3); % find analytic surface matrix
        [S2,R2,T2] = GET_SURFACE_TRANSFORMATIONS_SIM(cell2,.3);
        M1 = T1 * R1 * S1; % surface matrices
        M2 = T2 * R2 * S2;
        
        % Third transform the approximate surface of E2 into the analytic
        % coords of E1
        r2_E1 = M1 \ r2;        % all the points defining E2, in E1 coords
        D = vecnorm(r2_E1(1:3,:)); % distance from origin in E1 space
        ix = find(D < 1); % any points less than 1 are within the surface of E1
        
        % Main calculation:
        if length(ix) < 10 % they do not overlap
            overlap = []; % add nothing to the overlaps
            COM     = [];
            Fpt_1   = [];
            direction = [];
        elseif (10 <= length(ix)) && (length(ix) < 100) % they only overlap by a little bit
            r1_E2   = M2 \ r1;
            D_help  = vecnorm(r1_E2(1:3,:));
            ix_E2   = find(D_help < 1);
            inters  = [r1(1:3,ix_E2), r2(1:3,ix)]';
            [~,vol] = convhull(inters(:,1),inters(:,2),inters(:,3));
            overlap = [nbors(j); vol];
            COM     = mean(inters)';
            Fpt_1   = get_force_point(COM, cell1);% Find force point on the surface of the cell of interest
            direc   = M1 * [0;0;0;1];
            direction = direc(1:3) - Fpt_1;
            direction = direction/norm(direction);
        else % they overlap by a significant amount
            v1       = [x1,y1,z1]; % vertices defining cell1
            v2       = [x2,y2,z2]; % vertices defining cell2
			r1_E2    = M2 \ r1;
			D_help   = vecnorm(r1_E2(1:3,:));
            ix_E2    = find(D_help < 1);
            inters   = [r1(1:3,ix_E2), r2(1:3,ix)]';
            COM      = mean(inters)'; % a point common to both cells
			[inters] = INTERSECTION(v1, v2, COM');
            [~,vol]  = convhull(inters(:,1), inters(:,2), inters(:,3)); % returns the volume of overlap
            overlap  = [nbors(j); vol]; % add to overlaps array
            COM      = mean(inters)'; % center of the overlap volume
            Fpt_1    = get_force_point(COM, cell_list(ix_o_i));
            direc    = M1 * [0;0;0;1];
            direction = direc(1:3) - Fpt_1;
            direction = direction/norm(direction);
            % Find force point on the surface of the cell of interest:
        end
        % Store info to cell array:
        Force_pts = [Force_pts, Fpt_1];
        Overlap_pts = [Overlap_pts, COM];
        Overlaps = [Overlaps, overlap];
        Directions = [Directions, direction];
%         cell_list(ix_o_i).Force_pts   = [cell_list(ix_o_i).Force_pts, Fpt_1];
%         cell_list(ix_o_i).Overlap_pts = [cell_list(ix_o_i).Overlap_pts, COM];
%         cell_list(ix_o_i).Overlaps    = [cell_list(ix_o_i).Overlaps, overlap];
    end
    
%     % Trim to unique overlap readings:
%     if ~isempty(cell_list(ix_o_i).Overlaps)
%         [~,ix_unique,~] = unique(cell_list(ix_o_i).Overlaps(1,:));
%         cell_list(ix_o_i).Overlaps    = cell_list(ix_o_i).Overlaps(:,ix_unique);
%         cell_list(ix_o_i).Overlap_pts = cell_list(ix_o_i).Overlap_pts(:,ix_unique);
%         cell_list(ix_o_i).Force_pts   = cell_list(ix_o_i).Force_pts(:,ix_unique);
%     end

end

function [m_pts, d_pts] = CHITIN_INTERACTIONS(cell_list, n)
% Obtain all the chitin interactions acting on cell n:

    % Mother interaction:
    mother_cell = cell_list(n).Mother;
    if ~isempty(mother_cell) % if this cell has a mother
        cell_1 = cell_list(mother_cell + 1);
        cell_2 = cell_list(n);
        [Sm, Rm, Tm] = GET_SURFACE_MATRICES(cell_1, 0);
        [Sd, Rd, Td] = GET_SURFACE_MATRICES(cell_2, 0);
        ix = find(cell_1.Daughters == n - 1);
        pt_1 = Tm * Rm * [cell_1.BudXYZ(:,ix); 1]; % point on mother cell surface
        pt_2 = Td * Rd * Sd * [-1; 0; 0; 1]; % point on daughter cell surface
        m_pts = pt_1(1:3) - pt_2(1:3); % displacement vector
    else
        m_pts = [0;0;0];
    end
    
    % Daughter interactions:
    daughter_cells = cell_list(n).Daughters;
    d_pts = [];
    if ~isempty(daughter_cells) % if this cell has any daughters
        for d = 1:length(daughter_cells)
            cell_1 = cell_list(n);
            cell_2 = cell_list(daughter_cells(d)+1);
            [Sm, Rm, Tm] = GET_SURFACE_MATRICES(cell_1, 0);
            [Sd, Rd, Td] = GET_SURFACE_MATRICES(cell_2, 0);
            pt_1 = Tm * Rm * [cell_1.BudXYZ(:,d); 1];
            pt_2 = Td * Rd * Sd * [-1; 0; 0; 1];
            d_pts = [d_pts, pt_2(1:3) - pt_1(1:3)];
        end
    else
        d_pts = [0;0;0];
    end
end

function [Fpt] = get_force_point(COM, cell_o_i)
    % Given a point that lies within the ellipsoid, extrapolate out to the
    % surface of the ellipsoid.
    [S,R,T] = GET_SURFACE_MATRICES(cell_o_i,0);
    M       = T * R * S;
    COM_E   = M \ [COM; 1];
    Fpt_E   = [COM_E(1:3)/norm(COM_E(1:3)); 1];
    Fpt_g   = M * Fpt_E;
    Fpt     = Fpt_g(1:3);

end