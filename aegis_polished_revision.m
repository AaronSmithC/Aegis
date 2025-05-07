% Project Aegis - Clean Version
clf
hold on
axis equal
axis off

%% 1. Define parameters
R1 = 200; % Major radius of torus
R2 = 50;  % Minor radius (tube radius)
R1_small = 200; % Major radius for smaller torus
R2_small = 40;  % Minor radius for smaller torus
R1_shield = 200; % Major radius of the water shield boundary
R2_shield = 47; %three meter shield width

u = linspace(0, 2*pi, 50);
v = linspace(0, 2*pi, 50);
[u, v] = meshgrid(u, v);

% Define three ring levels (z-separated)
Z_offset = [0, 200, 400];

%% 2. Create 3 Large Rings
for i = 1:3
    X = (R1 + R2*cos(v)) .* cos(u);
    Y = (R1 + R2*cos(v)) .* sin(u);
    Z = R2*sin(v) + Z_offset(i);
    surf(Z, X, Y, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
end

%% 3. Create 3 Smaller Rings
for i = 1:3
    x_small = (R1_small + R2_small*cos(v)) .* cos(u);
    y_small = (R1_small + R2_small*cos(v)) .* sin(u);
    z_small = R2_small*sin(v) + Z_offset(i);
    surf(z_small, x_small, y_small, 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'none');
end

%Water Shield Boundary
%% 2. Create water shields
for i = 1:3
    X = (R1_shield + R2_shield*cos(v)) .* cos(u);
    Y = (R1_shield + R2_shield*cos(v)) .* sin(u);
    Z = R2_shield*sin(v) + Z_offset(i);
    surf(Z, X, Y, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
end


%% 4. Central Axis (cylinder)
rc = 20;   % radius of cylinder
hc = 650;  % height
[Xc, Yc, Zc] = cylinder(rc);
Zc = Zc * hc - 125;
surf(Zc, Xc, Yc, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

%% 5. Spokes
num_spokes = 5;
spoke_length = 150;
for i = 1:num_spokes
    angle = (i-1) * 2*pi / num_spokes;
    x = spoke_length * cos(angle);
    y = spoke_length * sin(angle);

    for j = 0:2 % three rings
        line([Z_offset(j+1), Z_offset(j+1)], [0, x], [0, y], 'Color', 'black', 'LineWidth', 5);
    end
end

%% 1. First, create the Docking Ring (same as before)
ringRadius = 80;  % Radius of docking ring
ringTubeRadius = 10; % Thickness of ring

theta = linspace(0, 2*pi, 50);
phi = linspace(0, 2*pi, 30);
[theta, phi] = meshgrid(theta, phi);

X_ring = (ringRadius + ringTubeRadius*cos(phi)) .* cos(theta);
Y_ring = (ringRadius + ringTubeRadius*cos(phi)) .* sin(theta);
Z_ring = ringTubeRadius * sin(phi);

centerAxisX = 0;
centerAxisY = 0;
centerAxisZ = 800; % height of station center
dockOffset = 350; % distance along Z

surf(Z_ring + centerAxisZ - dockOffset, ...
     X_ring + centerAxisX, ...
     Y_ring + centerAxisY, ...
     'FaceColor', [0.5 0.5 0.7], 'EdgeColor', 'none');
hold on

%% 1. Docking Ring (same)
ringRadius = 80;  
ringTubeRadius = 10; 

theta = linspace(0, 2*pi, 50);
phi = linspace(0, 2*pi, 30);
[theta, phi] = meshgrid(theta, phi);

X_ring = (ringRadius + ringTubeRadius*cos(phi)) .* cos(theta);
Y_ring = (ringRadius + ringTubeRadius*cos(phi)) .* sin(theta);
Z_ring = ringTubeRadius * sin(phi);

centerAxisX = 0;
centerAxisY = 0;
centerAxisZ = 800; 
dockOffset = 350; 

surf(Z_ring + centerAxisZ - dockOffset, ...
     X_ring + centerAxisX, ...
     Y_ring + centerAxisY, ...
     'FaceColor', [0.5 0.5 0.7], 'EdgeColor', 'none');
hold on

%% 2. Fancier Ships!

% Base ship body dimensions
bodyLength = 20;
bodyRadius = 8;

% Nose cone dimensions
noseLength = 10;  % shorter sharp nose

% Fin dimensions
finSpan = 6; 
finLength = 8;

% Ship location setup
numShips = 4;
for i = 1:numShips
    angle = (i-1) * 2*pi / numShips; % 0, 90, 180, 270 degrees

    % Ship center
    shipPosX = centerAxisX + (ringRadius + 20) * cos(angle);
    shipPosY = centerAxisY + (ringRadius + 20) * sin(angle);
    shipPosZ = centerAxisZ - dockOffset;

    %% Body - cylinder
    [shipX, shipY, shipZ] = cylinder(bodyRadius, 20);
    shipZ = shipZ * bodyLength;

    shiftedBodyX = shipX + shipPosX;
    shiftedBodyY = shipY + shipPosY;
    shiftedBodyZ = shipZ + shipPosZ;

    if abs(cos(angle)) > 0.5
        surf(shiftedBodyZ, shiftedBodyX, shiftedBodyY, ...
             'FaceColor', [0.6 0.6 1], 'EdgeColor', 'none');
    else
        surf(shiftedBodyZ, shiftedBodyY, shiftedBodyX, ...
             'FaceColor', [0.6 0.6 1], 'EdgeColor', 'none');
    end

    %% Nose - cone
    [noseX, noseY, noseZ] = cylinder([bodyRadius 0], 20);
    noseZ = noseZ * noseLength;

    shiftedNoseX = noseX + shipPosX;
    shiftedNoseY = noseY + shipPosY;
    shiftedNoseZ = noseZ + shipPosZ + 20; % offset so tip points center

    if abs(cos(angle)) > 0.5
        surf(shiftedNoseZ, shiftedNoseX, shiftedNoseY, ...
             'FaceColor', [0.9 0.9 1], 'EdgeColor', 'none');
    else
        surf(shiftedNoseZ, shiftedNoseY, shiftedNoseX, ...
             'FaceColor', [0.9 0.9 1], 'EdgeColor', 'none');
    end

    %{
%% Fins (Corrected)
% 4 simple triangular fins
finX = [0 0 -finLength];
finY = [0 finSpan finSpan/2];
finZ = [0 0 0]; % flat fins along Z=0

% Repeat fins rotated around body
finAngles = [0 pi/2 pi 3*pi/2];
for j = 1:4
    rot = finAngles(j);

    % Rotate fin shape around the body
    R = [cos(rot) -sin(rot); sin(rot) cos(rot)];
    newFin = R * [finX; finY];

    % Build the fin's vertices
    finXr = newFin(1,:) + shipPosX;
    finYr = newFin(2,:) + shipPosY;
    finZr = shipPosZ + finZ; % match size

    % Now draw it correctly
    fill3(finZr, finXr, finYr, [0.7 0.7 0.7], 'EdgeColor', 'none');
end
    %}
end

hold on;


%% 7. Rivets (small spheres)
numRivetsTheta = 20;
numRivetsPhi = 10;
thetaR = linspace(0, 2*pi, numRivetsTheta);
phiR = linspace(0, 2*pi, numRivetsPhi);
[thetaR, phiR] = meshgrid(thetaR, phiR);

xRivet = (R1 + R2*cos(phiR)) .* cos(thetaR);
yRivet = (R1 + R2*cos(phiR)) .* sin(thetaR);
zRivet = R2 * sin(phiR);

rivetSize = 3;
for i = 1:numel(xRivet)
    [sx, sy, sz] = sphere(10);
    for j = 0:2 % three rings
        surf(rivetSize*sx + zRivet(i) + Z_offset(j+1), ...
             rivetSize*sy + xRivet(i), ...
             rivetSize*sz + yRivet(i), ...
             'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    end
end

%% Add Base Disc for Communication Arrays

% Parameters for base disc
baseRadius = 50;   % Bigger than the dish attachment points
baseThickness = 5; % Thickness of the disc
baseColor = [0.5 0.5 0.7]; % Light metallic gray-blue

% Create base for top
[baseX, baseY, baseZ] = cylinder(baseRadius, 50);
baseZ = baseZ * baseThickness; % scale to thickness
baseZ = baseZ + (hc - 125 + 5); % shift up to top

surf(baseZ, baseX, baseY, ...
    'FaceColor', baseColor, 'EdgeColor', 'none');

% Create base for bottom
[baseX2, baseY2, baseZ2] = cylinder(baseRadius, 50);
baseZ2 = baseZ2 * baseThickness; % scale to thickness
baseZ2 = baseZ2 + (-125 - 5 - baseThickness); % shift to bottom

surf(baseZ2, baseX2, baseY2, ...
    'FaceColor', baseColor, 'EdgeColor', 'none');
%% Multiple Communications Dishes at End of Central Hub - Top
arrayBaseZ_top = hc - 125 + 40; % Top end
arrayBaseX = 0;
arrayBaseY = 0;

% Dish parameters
dishRadius = 30;
dishDepth = 10;
numDishes = 4;
dishAngleStep = 2*pi/numDishes;

for i = 1:numDishes
    angle = (i-1) * dishAngleStep;
    offsetDistance = 40;
    centerX = arrayBaseX + offsetDistance * cos(angle);
    centerY = arrayBaseY + offsetDistance * sin(angle);
    centerZ = arrayBaseZ_top;
    
    [theta, rho] = meshgrid(linspace(0, 2*pi, 30), linspace(0, 1, 10));
    Xdish = rho .* dishRadius .* cos(theta);
    Ydish = rho .* dishRadius .* sin(theta);
    Zdish = -(rho.^2) * dishDepth;
    
    tiltAngle = 25; 
    tiltRad = deg2rad(tiltAngle);
    R = [cos(angle) -sin(angle) 0;
         sin(angle)  cos(angle) 0;
         0           0          1];
    T = [1 0 0;
         0 cos(tiltRad) -sin(tiltRad);
         0 sin(tiltRad)  cos(tiltRad)];
    
    points = [Xdish(:), Ydish(:), Zdish(:)] * T' * R';
    
    XdishR = reshape(points(:,1), size(Xdish)) + centerX;
    YdishR = reshape(points(:,2), size(Ydish)) + centerY;
    ZdishR = reshape(points(:,3), size(Zdish)) + centerZ;
    
    surf(ZdishR, XdishR, YdishR, ...
        'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');
    
    rodHeight = 20;
    rodRadius = 1.5;
    [rodX, rodY, rodZ] = cylinder(rodRadius, 20);
    rodZ = rodZ * rodHeight;
    rodX = rodX + centerX;
    rodY = rodY + centerY;
    rodZ = rodZ + centerZ;
    surf(rodZ, rodX, rodY, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
end

%% Now, Communications Dishes at Bottom End
arrayBaseZ_bottom = -125 - 40; % Bottom end

for i = 1:numDishes
    angle = (i-1) * dishAngleStep;
    offsetDistance = 40;
    centerX = arrayBaseX + offsetDistance * cos(angle);
    centerY = arrayBaseY + offsetDistance * sin(angle);
    centerZ = arrayBaseZ_bottom;
    
    [theta, rho] = meshgrid(linspace(0, 2*pi, 30), linspace(0, 1, 10));
    Xdish = rho .* dishRadius .* cos(theta);
    Ydish = rho .* dishRadius .* sin(theta);
    Zdish = -(rho.^2) * dishDepth;
    
    tiltAngle = -25; % opposite tilt
    tiltRad = deg2rad(tiltAngle);
    R = [cos(angle) -sin(angle) 0;
         sin(angle)  cos(angle) 0;
         0           0          1];
    T = [1 0 0;
         0 cos(tiltRad) -sin(tiltRad);
         0 sin(tiltRad)  cos(tiltRad)];
    
    points = [Xdish(:), Ydish(:), Zdish(:)] * T' * R';
    
    XdishR = reshape(points(:,1), size(Xdish)) + centerX;
    YdishR = reshape(points(:,2), size(Ydish)) + centerY;
    ZdishR = reshape(points(:,3), size(Zdish)) + centerZ;
    
    surf(ZdishR, XdishR, YdishR, ...
        'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');
    
    rodHeight = 20;
    rodRadius = 1.5;
    [rodX, rodY, rodZ] = cylinder(rodRadius, 20);
    rodZ = rodZ * rodHeight;
    rodX = rodX + centerX;
    rodY = rodY + centerY;
    rodZ = rodZ + centerZ;
    surf(rodZ, rodX, rodY, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
end


%% 8. Final touches
grid off;
shading interp;
camlight('headlight')
lighting gouraud
material shiny
set(gca, 'Color', 'k')
view(45, 30)
title('Project Aegis - Simple Model','Color','w','FontSize',16)

%% Set up camera view and lighting for nice render
axis equal
grid off
shading interp
camlight('headlight')
lighting gouraud
material shiny
set(gca, 'Color', 'k') % Black background
view(45, 30) % Good 3D angle
title('Aegis Station','Color','black','FontSize',18)
