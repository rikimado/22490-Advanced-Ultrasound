function plotTimeDomain(t_exc, excitation, t_imp, impulse_response, fs, f0)
% plotTimeDomain - Plots excitation, impulse response, and windowed impulse in time domain
%
% Inputs:
%   excitation       - vector containing the excitation waveform
%   impulse_response - vector containing the impulse response waveform
%   fs               - sampling frequency (Hz)
%   f0               - center frequency (Hz) for marking cycles

% Compute windowed impulse response
impulse_windowed = impulse_response .* hann(numel(impulse_response))';

% Plot signals
figure;
plot(t_exc*1e6, excitation, 'k', 'LineWidth', 1.2, 'DisplayName', 'e(t) (1 cycle) ');
hold on;
plot(t_imp*1e6, impulse_response, 'r--', 'LineWidth', 1.2, 'DisplayName', 'v(t) (2 cycles)');
plot(t_imp*1e6, impulse_windowed, 'g', 'LineWidth', 1.2, 'DisplayName', 'windowed v(t)');

% Number of cycles (approx) for markers
N_exc = round(length(excitation) / (fs/f0)); 
N_imp = round(length(impulse_response) / (fs/f0));

% Add vertical dashed lines at end of each cycle
for i = 1:N_exc
    xline(i*(1/f0)*1e6, 'b--', 'HandleVisibility','off');  % Excitation cycles
end
for i = 1:N_imp
    xline(i*(1/f0)*1e6, 'm--', 'HandleVisibility','off');  % Impulse cycles
end

xlabel('Time (µs)');
ylabel('Amplitude');
legend();
title('Time Domain Waveforms');
grid on;

end