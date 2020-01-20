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
% Proofs if the wavespeeds for a specific source was used before.
% version='$Rev: 8 $ ($Date: 2017-05-29 13:51:14 +0200 (Mo, 29 Mai 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

source_frequency_vec = source_matrix(:,1);
source_type_vec = source_matrix(:,2);
source_matrix_tmp = [source_frequency_vec(1:source-1), source_type_vec(1:source-1)];
used_source = source_matrix(source,:);
check = ismember(used_source,source_matrix_tmp);

if check(1) == true && check(2) == true
    fprintf('The same source was used before. \n');
    checking_source = true;
end
