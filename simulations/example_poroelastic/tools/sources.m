%-----------------------------------------------------------------------
%   Copyright 2015-2017 Janis Heuel (Ruhr-Universität Bochum, GER)
%   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
%
%   This file is part of NEXD 1D.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with NEXD 1D. If not, see <http://www.gnu.org/licenses/>.
%---------------------------------
% % version='$Rev: 3 $ ($Date: 2016-08-08 14:28:52 +0100 (Mo, 08 August 2016) $, $Author: Janis Heuel, Marc S. Boxberg $)'
% Thsi script reads the parameters from 'sources' and calculates the
% the frequency range of three wavelet types:
% Ricker:   1
% Green:    2
% External: 3
%

cd (workpath)
cd ../input
figure_number = figure_number+source;


%% Reading the parameters from 'sources' for the different wavelets
read_sources                 = table2cell(readtable('sources'));
find_BEGIN                   = find(strcmp(read_sources, 'BEGIN'));
find_END                     = find(strcmp(read_sources, 'END'));

del_lines_in_sources         = read_sources(find_BEGIN+source:find_END-1);
convert_sources              = cell2mat(del_lines_in_sources(1));
ricker_wavelet               = strfind(convert_sources, 'ricker');   % Ricker wavelet
green_wavelet                = strfind(convert_sources, 'green');    % Green wavelet
external_wavelet             = strfind(convert_sources, 'external'); % External wavelet
del_end_of_convert_sources   = convert_sources(1:end-1);

if ricker_wavelet == 1 % Number of del_begin depends on length of word (e.g. ricker has 6 letters)
    del_begin = 6;
elseif green_wavelet == 1
    del_begin = 5;
elseif external_wavelet == 1
    del_begin = 8;
else
    fprintf('Unknown wavelet. The script stops here. \n');
    fprintf('Possible types of wavelet: ricker, green, external \n');
    cd (workpath)
    cd ../tools
    return
end

del_begin_of_convert_sources = del_end_of_convert_sources(1+del_begin:end);
convert_del_begin_to_array   = str2num(del_begin_of_convert_sources);
f_c                          = convert_del_begin_to_array(2);  % Most energetic frequency (in Hertz) for the wavelet
omega_p                      = 2*pi*f_c;                       % Most energetic angular frequency (in radians per second)


if ricker_wavelet == 1
    source_type = 1;

    % Fourier transforamtion of the Ricker wavelet
    omega_ricker   = linspace(0,25*f_c,1e6);
    Fourier_ricker = (2.*omega_ricker.^2) ./ (sqrt(pi)*omega_p^3) .* (exp(-(omega_ricker.^2) ./(omega_p^2))); % analytical fourier transformation of the Ricker wavelet

    % Finding roots of F for maximal value of the angular frequency
    roots_Fourier_ricker                      = find(Fourier_ricker <= wavelet_amplitude*max(Fourier_ricker));
    derivative_of_roots_Fourier_ricker        = diff(roots_Fourier_ricker);
    find_discontinuty                         = find(derivative_of_roots_Fourier_ricker==max(derivative_of_roots_Fourier_ricker))+1;
    find_entry_in_roots_Fourier_ricker        = roots_Fourier_ricker(find_discontinuty);
    omega_max                                 = omega_ricker(find_entry_in_roots_Fourier_ricker); % Determination of maximal value of the angular frequency

elseif green_wavelet == 1
    source_type = 2;
    % Fourier transformation of the green wavelet
    omega_green = linspace(0,25*f_c,1e6);
    Fourier_green = (exp(-omega_green.^2 ./ (4*pi^2*f_c^2))) / (sqrt(2*pi)); % analytical fourier transformation of green wavelet

    % Finding roots of the green fourier transform
    roots_Fourier_green = find(Fourier_green <= wavelet_amplitude*max(Fourier_green));
    omega_max = omega_green(roots_Fourier_green(1));

elseif external_wavelet == 1
    source_type = 3;
    external_table = str2num(char(table2array(readtable(char(extwavelet))))); % Reading the .txt file which includes the external wavelet
    time           = external_table(:,1);                                     % Time vector of the external wavelet
    amplitude      = external_table(:,2);                                     % Amplitude vector of the external wavelet

    sampling_rate  = length(amplitude)/(max(time)-min(time)); % sampling rate calculated by the number of samples and the length of the external wavelet

    trans_length  = pow2(nextpow2(length(amplitude)));                % Transformation of the length
    fast_ft       = fft(amplitude,trans_length);                      % Fast Fourier Transform
    frequency     = (0:trans_length-1)*(sampling_rate/trans_length);  % Frequency range
    power         = fast_ft.*conj(fast_ft)/trans_length;              % Power of the fourier transform

    % Parameters for roots of external wavelet
    angular_frequency_roots = 2*pi.*frequency(1:floor(trans_length/2)); % Frequency in rad/s
    power_roots             = power(1:floor(trans_length/2));           % power of the external wavelet

    % Finding roots of F for maximal value of the angular frequency
    roots_Fourier_external                      = find(power_roots <= wavelet_amplitude*max(power_roots));
    derivative_of_roots_Fourier_external        = diff(roots_Fourier_external);
    if isempty(find(derivative_of_roots_Fourier_external>1, 1))
        find_discontinuty = 1;
    else
        find_discontinuty                       = find(derivative_of_roots_Fourier_external==max(derivative_of_roots_Fourier_external))+1;
    end
    find_entry_in_roots_Fourier_external        = roots_Fourier_external(find_discontinuty);
    omega_max                                   = angular_frequency_roots(find_entry_in_roots_Fourier_external); % Determination of maximal value of the angular frequency
end

source_matrix(source,:) = [f_c source_type]; % Matrix with entries of each source


%% Printing wavelet type
fprintf('\n');
fprintf('--------------------------------------------------- \n');
if ricker_wavelet == 1
    fprintf('Type of the wavelet: Ricker \n');
elseif green_wavelet == 1
    fprintf('Type of the wavelet: Green \n')
elseif external_wavelet == 1
    fprintf('Type of the wavelet: External \n');
end
fprintf('Most energetic frequency of the wavelet: %g Hz \n', f_c);
fprintf('Frequency spectrum: 0 - %g rad/s \n', omega_max);
fprintf('Frequency spectrum: 0 - %g Hz \n', omega_max/(2*pi));
fprintf('--------------------------------------------------- \n');
fprintf('\n');
