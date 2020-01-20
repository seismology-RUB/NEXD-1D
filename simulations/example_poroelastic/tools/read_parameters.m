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
% This script reads the parameters from 'mat_bap' and parfile.
% For reading the parameters from the files, the function 'read_input.m'
% is used. For more information see description of this dunction.
% version='$Rev: 8 $ ($Date: 2017-05-29 13:51:14 +0200 (Mo, 29 Mai 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

%% Checking if calculation is elastic or poroelastic.
poroelastic = read_input(str2func('../input/parfile'),str2func('poroelastic'),1,0);
fprintf('\n')
if poroelastic == 1
    fprintf('Calculation is poroelastic. \n')
elseif poroelastic == 0
    fprintf('Calcualtion is elastic. \n')
else
    fprintf('Error in parfile at point poroelastic!\n');
    cd (workpath)
    cd ../tools
    return
end

%% Reading the number of materials
number_of_materials = read_input(str2func('../input/parfile'),str2func('matn'),0,0);
if number_of_materials <= 0
    fprintf('Number of materials is zero. The script stops here! \n');
    cd (workpath)
    cd ../tools
    return
elseif number_of_materials > 0
    fprintf('Number of materials: %g \n', number_of_materials);
end


%% Reading the number of fluids
number_of_fluids = read_input(str2func('../input/parfile'),str2func('fluidn'),0,0,0);
if ~ ismember(number_of_fluids,[0,1,2])
    fprintf('Number of fluids is greater than 2. The script stops here. \n')
    cd (workpath)
    cd ../tools
    return
end
number_of_fluids_save = number_of_fluids;
fprintf('\n');
fprintf('Number of fluids: %g \n',number_of_fluids);


%% Reading xmin
xmin = read_input(str2func('../input/parfile'),str2func('xmin'),0,0,0);
fprintf('Minimal x-value: %g \n', xmin);


%% Reading number of sources
srcn = read_input(str2func('../input/parfile'),str2func('srcn'),0,0,0);
fprintf('Number of sources: %g \n', srcn);
% Allocate field for sources
source_vector = zeros(srcn,2);  % -> is used in proof_source.m


% Finding external wavelet in 'parfile'
extwavelet = read_input(str2func('../input/parfile'),str2func('extwavelet'),0,0,0);


%% Checking if tortuosity is given in 'matpropporo' or if tortuosity
%% will be calculated by a further parameter r.
tortuosity = read_input(str2func('../input/parfile'),str2func('calculate_tortuosity'),1,0,0);
if tortuosity == 1
    calculate_tortuosity = true;
elseif tortuosity == 0
    calculate_tortuosity = false;
else
    fprintf('calculate_tortuosity has to be set as true or false in parfile. \n');
    cd (workpath)
    cd ../tools
    return
end


cd (workpath)
cd ../input


%% Reading materials from 'mat_bap'
read_mat_bap      = table2cell(readtable('mat_bap'));
find_BEGIN_porous = strfind(read_mat_bap,'BEGIN');
b_loop = 0;
l_loop = 1;
while b_loop ~= 1 % finds the line with the keyword 'BEGIN'
    a_loop  = cell2mat(find_BEGIN_porous(l_loop));
    b_loop  = length(a_loop);
    l_loop  = l_loop+1;
end
find_BEGIN_porous = l_loop-1;

find_END_porous = strfind(read_mat_bap,'END');
b_loop = 0;
l_loop = 1;
while b_loop ~= 1 % finds the line with the keyword 'END'
    a_loop  = cell2mat(find_END_porous(l_loop));
    b_loop  = length(a_loop);
    l_loop  = l_loop+1;
end
find_END_porous = l_loop-1;
mat_boundaries  = zeros(find_END_porous-(find_BEGIN_porous+1),1); % Boundaries of the different materials
mat_types       = zeros(find_END_porous-(find_BEGIN_porous+1),1); % Order of material types
for i = 1:(find_END_porous-(find_BEGIN_porous+1))
    mat_bap_cell      = cell2mat(read_mat_bap(find_BEGIN_porous+i));
    mat_bap_vec       = str2num(mat_bap_cell);
    mat_boundaries(i) = mat_bap_vec(1); % Generating the vectors for materials
    mat_types(i)      = mat_bap_vec(2);
end
cd (workpath)


