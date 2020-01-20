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
%% This cript plots the output of calculation of wavespeed and inverse quality factor.
% version='$Rev: 8 $ ($Date: 2017-05-29 13:51:14 +0200 (Mo, 29 Mai 2017) $, $Author: Janis Heuel, Marc S. Boxberg $)'

%##########################################################
% YOU CAN CHANGE FOR OTHER AXIS LABELS AND UNITS THE FOLLOWING LINES
xlabel_unit = Frequency;          % unit of x-axis
xaxis_label = 'Frequency f [Hz]'; % label of x-axis

%xlabel_unit = omega;
%xaxis_label = 'Angular frequency \omega [rad/s]';

y_axis_label_wavespeed = 'Wavespeed [m/s]';
y_axis_label_Q         = 'Q^{-1}';
%##########################################################

figure(figure_number)
set(figure(figure_number),'name',['Wavespeeds of material ' num2str(u) ' and source ' num2str(source)])
if attenuation == true
    figure(figure_number+1)
    set(figure(figure_number+1),'name',['Attenuation of material ' num2str(u) ' and source ' num2str(source)])
end

if number_of_fluids == 0 || phi_in == 0
    figure(figure_number)
        semilogx(xlabel_unit,v1)
        legend('P')
        ylabel(y_axis_label_wavespeed)
        xlabel(xaxis_label)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(source)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            semilogx(xlabel_unit,Q1)
            legend('P')
            if save_plot == true
                print(['a' num2str(u) 's' num2str(source)],formattype)
            end
        end
elseif number_of_fluids == 1
    for i = 1:2
        figure(figure_number)
        subplot(2,1,i)
        if i == 1
            semilogx(xlabel_unit,v1)
            legend('P1')
        elseif i == 2
            semilogx(xlabel_unit,v2)
            legend('P2')
        end
        ylabel(y_axis_label_wavespeed)
        xlabel(xaxis_label)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(source)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            subplot(2,1,i)
            if i == 1
                semilogx(xlabel_unit,Q1)
                legend('P1')
            elseif i == 2
                semilogx(xlabel_unit,Q2)
                legend('P2')
            end
            ylabel(y_axis_label_Q)
            xlabel(xaxis_label)
            if save_plot == true
                print(['a' num2str(u) 's' num2str(source)],formattype)
            end
        end
    end
elseif number_of_fluids == 2
    for i = 1:3
        figure(figure_number)
        subplot(3,1,i)
        if i == 1
            semilogx(xlabel_unit,v1)
            legend('P1')
        elseif i == 2
            semilogx(xlabel_unit,v2)
            legend('P2')
        else
            semilogx(xlabel_unit,v3)
            legend('P3')
        end
        ylabel(y_axis_label_wavespeed)
        xlabel(xaxis_label)
        if save_plot == true
            print(['w' num2str(u) 's' num2str(source)],formattype)
        end
        if attenuation == true
            figure(figure_number+1)
            subplot(3,1,i)
            if i == 1
                semilogx(xlabel_unit,Q1)
                legend('P1')
            elseif i == 2
                semilogx(xlabel_unit,Q2)
                legend('P2')
            elseif i == 3
                semilogx(xlabel_unit,Q3)
                legend('P3')
            end
            ylabel(y_axis_label_Q)
            xlabel(xaxis_label)
            if save_plot == true
                print(['a' num2str(u) 's' num2str(source)],formattype)
            end
        end
    end
end
