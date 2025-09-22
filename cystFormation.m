% Create a computer model of a cyst phantom. The phantom contains
% N point targets separated by 5 mm and a 10 mm water filled cyst.
% All scatterers are situated in a box of (x,y,z) mm.
%
% Calling: [positions, amp] = cystFormation(N);
%
% Input: N - Number of scatterers in the phantom
%
% Output:
% - positions: positions of the scatterers.
% - amp: amplitude of the scatterers.
%

%%

function [positions, amp] = cystFormation(N)
x_size = 15/1000;    % Width of box [m]
y_size = 0/1000;     % Transverse width of box [m]
z_size = 15/1000;    % Height of box [m]
z_start = 22.5/1000;    % Start of box surface [m];

% Create the general scatterers:
% Generates a vector of N random numbers from a uniform 
% distribution between 0 and 1. Subtracting 0.5 
% from each element shifts the range from [0, 1) to [-0.5, 0.5).
% Multiplying the shifted values by a size scales the range to that size
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;

% Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);

% Make the water-filled cyst and set the amplitudes to zero inside
r=3/1000;   % Radius of cyst [m]
xc=0/1000;  % center of cyst [m]
zc=30/1000;
inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
amp = amp .* (1-inside);

% Create scatterers at fixed points for testing.
% dz=z_size/10;
% for i=N-9:N
%     x(i) = -15/1000;
%     y(i) = 0;
%     z(i) = z_start + (i-N+9)*dz;
%     amp(i) = 100;
% end
% 
% % Return the variables
positions=[x y z];
end