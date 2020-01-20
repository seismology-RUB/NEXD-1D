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
%% This script plots the output of calculation of wavespeed and inverse quality factor.
% version='$Rev: 8 $ ($Date: 2017-05-29 13:51:14 +0200 (Mo, 29 Mai 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

fprintf('\n');
fprintf('--------------------------------------------------- \n');
fprintf('Boundary of material %g: %g \n',u,mat_boundaries(u));
if number_of_fluids == 2
    if ismember(used_model,[1,3,5]) == true
        fprintf('Used permeability model: van Genuchten (1980) \n');
    elseif ismember(used_model,[2,4,6]) == true
        fprintf('Used permeability model: Brooks & Corey (1964) \n');
    end
end
if number_of_fluids == 0
    fprintf('Wavespeed of P for material %g: %g m/s \n', u, max(v1));
elseif number_of_fluids > 0
    fprintf('--------------------------------------------------- \n');
    fprintf('Maximal wavespeed of P1 for material %g: %g m/s \n', u, max(v1));
    fprintf('Minimal wavespeed of P1 for material %g: %g m/s \n', u, min(v1));
    fprintf('Wavespeed at the most energetic frequency of P1: %g m/s \n', v1_f_c(2));
    fprintf('--------------------------------------------------- \n');
    if attenuation == true
        fprintf('Maximal inverse quality factor of P1 for material %g: %.15g \n', u, max(Q1));
        fprintf('Minimal inverse quality factor of P1 for material %g: %.15g \n', u, min(Q1));
        fprintf('Inverse quality factor at the most energetic frequency of P1: %.15g \n', Q1_f_c(2));
        fprintf('--------------------------------------------------- \n');
    end
    if number_of_fluids == 1 || number_of_fluids == 2
        fprintf('Maximal wavespeed of P2 for material %g: %g m/s \n', u, max(v2));
        fprintf('Minimal wavespeed of P2 for material %g: %g m/s \n', u, min(v2));
        fprintf('Wavespeed at the most energetic frequency of P2: %g m/s \n', v2_f_c(2));
        fprintf('--------------------------------------------------- \n');
        if attenuation == true
            fprintf('Maximal inverse quality factor of P2 for material %g: %.15g \n', u, max(Q2));
            fprintf('Minimal inverse quality factor of P2 for material %g: %.15g \n', u, min(Q2));
            fprintf('Inverse quality factor at the most energetic frequency of P2: %.15g \n', Q2_f_c(2));
            fprintf('--------------------------------------------------- \n');
        end
    end
    if number_of_fluids == 2
        fprintf('Maximal wavespeed of P3 for material %g: %g m/s \n', u, max(v3));
        fprintf('Minimal wavespeed of P3 for material %g: %g m/s \n', u, min(v3));
        fprintf('Wavespeed at the most energetic frequency of P3: %g m/s \n', v3_f_c(2));
        fprintf('--------------------------------------------------- \n');
        if attenuation == true
            fprintf('Maximal inverse quality factor of P3 for material %g: %.15g \n', u, max(Q3));
            fprintf('Minimal inverse quality factor of P3 for material %g: %.15g \n', u, min(Q3));
            fprintf('Inverse quality factor at the most energetic frequency of P2: %.15g \n', Q3_f_c(2));
            fprintf('--------------------------------------------------- \n');
        end
    end
end
fprintf('\n');
