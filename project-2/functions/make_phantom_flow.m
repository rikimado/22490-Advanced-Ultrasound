function [x,y,z,amp,velocity] = make_flow_phantom(x_range, z_range, peak_velocity, num_scatterers, flowProfile)

% Generates random 3D positions (x,y,z) and corresponding flow velocites
% for a number of scatterers inside a vessel-shaped volume

% INPUTS:
% - x_range = range for scatterer distribution along the x-axis [m]
% - z_range = range for scatterer distribution along the z-axis [m]
% - num_scatterers = number of scatterers to generate.
% - peak_velocity = maximum flow velocity [m/s] used in the velocity profile.
% - flowProfle = a string (plug/parabolic), defining the flow velocity distribution.

% Outputs
% - x, y, z: 3D coordinates of the scatterers.
% - amp: random amplitude (scattering strength) of each scatterer (Gaussian-distributed).
% - velocity: flow velocity assigned to each scatterer, depending on the profile.

% scatterer position generation
x = x_range*(rand(1, num_scatterers)-0.5);
z = z_range*(rand(1, num_scatterers)-0.5);
y = zeros(1, num_scatterers);
amp = randn(1, num_scatterers);

switch (flowProfile)
    case 'plug'
        velocity = peak_velocity;
    case 'parabolic'
        r = z/z_range*2;
        velocity = peak_velocity*(1-r.^2);
end
