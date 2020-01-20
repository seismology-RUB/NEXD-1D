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
% Tool to calculate wavespeeds for porousmaterials in one dimension
% version='$Rev: 11 $ ($Date: 2018-07-17 10:49:11 +0200 (Di, 17 Jul 2018) $, $Author: Janis Heuel, Marc S. Boxberg $)'

clear all;
warning('off','all')
workpath = pwd;
cd (workpath)
cd ../
mkdir output
cd output/

%%###############################
%% CONVENIENCE PARAMETERS
plot_diag           = true;       % if plot_diag = true, the plots are displayed
save_plot           = false;      % if saves_plot = true, the plot is saved in available formattypes in matlab
formattype          = '-dpng';
attenuation_default = true;       % if attenuation_default = true, the invers quality factor is calculated
logfile             = true;       % Set true if logfile should be created, otherwise false
k_default           = 40;         % Length of all vectors
wavelet_amplitude   = 0.001;      % minimal value of amplitude of fourier transformation of the wavelet,
                                  % used for frequency estimation. E.g. wavelet_amplitude = 0.001, then 0.1 % of the
                                  % highest amplitude for the highest frequency
%%###############################
%% Logfile
if logfile
    diary on
    delete('logfile_calc_vel')
    diary('logfile_calc_vel')
end
save_plot_diag = plot_diag;
cd (workpath)
run('read_parameters.m')

%% CALCULATE VELOCITIES FOR EACH SOURCE

figure_number = 0; % figure_number is important to create for each iteration a figure
for source = 1:srcn
    cd (workpath)
    cd ../tools
    run('sources.m')

    checking_source = false;
    if source >= 2
        cd (workpath)
        cd ../tools
        run('proof_source.m')
    end

    if checking_source == false
        n_cell    = zeros(1,number_of_materials); % Zero vector for the number of elements
        n_cell_P1 = zeros(1,number_of_materials);
        n_cell_P2 = zeros(1,number_of_materials);
        n_cell_P3 = zeros(1,number_of_materials);

        %% CALCULATE VELOCITIES FOR EACH MATERIAL
        for u = 1:number_of_materials
            cd (workpath)
            cd ../tools
            run('porousmaterial.m')

            if phi_in == 0
                number_of_fluids = 0;
            end

            % Calculation of further parameters for calcualting wavespeeds
            if number_of_fluids == 0 % 0 fluids
                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        if phi_in == 0
                            b     = zeros(1,k); % Biot coefficient
                            inv_N = zeros(1,k); % Biot modulus
                        elseif phi_in ~= 0
                            b     = 1 - K_d ./ K_s;
                            inv_N = b ./ K_s;
                        end
                    elseif ismember(used_model,[5,6]) == true
                        b = 1/2 - sqrt(1/4 - K_d.*inv_N);
                    end
                    M        = 1./inv_N;
                    lambda_u = K_d + M.*b.^2 - 2/3.*my; % first Lamé Parameter (undrained)
                end

            elseif number_of_fluids == 1   % 1 fluid
                if S1_in ~= 0
                    S1 = 1.0 * ones(1,k);
                elseif S1_in == 0 && S2_in ~= 0
                    S1   = S2_in*ones(1,k);
                    K1   = K2;
                    rho1 = rho2;
                    eta1 = eta2;
                end
                phi1   = S1.*phi;
                rho    = (1-phi).*rhos + (phi1.*rho1);

                if ismember(used_model,[1,3,5]) == true   % van Genuchten
                    m = 1-1./fitting_n;
                end

                K_f = (S1./K1).^(-1);                         % Bulk modulus fluid

                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        b     = 1 - K_d ./ K_s;                 % Biot coefficient
                        inv_N = (b-phi) ./ K_s;                 % inverse Biot modulus
                    elseif ismember(used_model,[5,6]) == true
                        b = (phi+1)/2 - sqrt(((phi+1).^2)/4 - phi - K_d.*inv_N);
                    end
                    M = (inv_N + phi./K_f).^(-1);
                    lambda_u = K_d + M.*b.^2 - 2/3.*my; % first Lamè Parameter (undrained)
                end

                M = (inv_N + phi./K1).^(-1);

                omegac1 = eta1.*inv_T.*phi1 ./ (kappa.*rho1);
                kappa1  = kappa;

                % Implementing correction function
                if all(omegac1) == 0
                    C1 = ones(1,k);
                else
                    if eta1_in == 0
                        C1 = ones(1,k);
                    elseif eta1_in ~= 0
                        C1 = sqrt(1 + 0.5i .* omega ./ omegac1);
                    end
                end

                rho_star     = (1-phi).*rhos + (1-inv_T) .* phi1.*rho1;
                gamma        = phi1 .* inv_T;
                Lambda_star  = kappa1 ./ (C1.*phi1.*eta1);
                Lambda1tilde = rho_star.*Lambda_star ./ (phi1.*(inv_T.^2-inv_T));
                Lambda2tilde = (inv_T./(rho1.*Lambda_star) + 1./Lambda1tilde).^(-1);

            elseif number_of_fluids == 2 % 2 fluids
                phi1 = S1.*phi;
                phi2 = S2.*phi;
                rho  = (1-phi).*rhos + (phi1.*rho1 + phi2.*rho2);

                % Calculate effective saturations
                S1eff = (S1-Sr1)./(1-Sr2-Sr1);
                S2eff = (S2-Sr2)./(1-Sr1-Sr2);
                
                if ismember(used_model,[1,3,5]) == true
                    m = 1-1./fitting_n;
                    G = - (rho1.*g)./((fitting_n-1).*fitting_chi) .* ((S1eff.^(-1./m)-1).^(1./fitting_n-1)) .* S1eff.^(-(1./m+1)); % dp_c/dS_1 van Genuchten
                elseif ismember(used_model,[2,4,6]) == true
                    G = - p_b ./ (lambda_BC * (Sr2-Sr1)) .* S1eff.^(-(1+lambda_BC)./lambda_BC); % dp_c/dS_1 Brooks & Corey
                elseif ismember(used_model,[7,8,9]) == true
                    G = - 2 .* A .* (1 ./ ((S1-Sr1).^3) + Sr2.^2 ./ ((1-S1).^3 .* (1-Sr1-Sr2).^2)); % dp_c/dS_1 Douglas Jr. et al.
                end

                K_f = (S1./K1 + S2./K2).^(-1);

                if ismember(used_model,[3,4,5,6]) == true
                    if ismember(used_model,[3,4]) == true
                        b        = 1 - K_d ./ K_s;
                        inv_N    = (b-phi) ./ K_s;
                    elseif ismember(used_model,[5,6]) == true
                        b = (phi+1)/2 - sqrt(((phi+1).^2)/4 - phi - K_d.*inv_N);
                    end
                    M        = (inv_N + phi1./K1 + phi2./K2 - ((S1.*S2.*phi) ./ (K1.*K2) + S1.*S2.*inv_N .* (S1 ./ K2 + S2 ./ K1)).*G).^(-1);
                    lambda_u = K_d + (1+S1.*S2.*(S1./K2+S2./K1).*G).*M.*b.^2 - 2/3.*my;
                end

                % Calcualtion of realtive permeability
                if ismember(used_model,[1,3,5]) == true
                    % van Genuchten
                    kappa1rel = sqrt(S1eff) .* (1 - (1-S1eff.^(1./m)).^m).^2;
                    kappa2rel = sqrt(1-S1eff) .* (1 - S1eff.^(1./m)).^(2.*m);
                elseif ismember(used_model,[2,4,6]) == true
                    % Brooks & Corey
                    kappa1rel = S1eff.^((2+3.*lambda_BC)./lambda_BC);
                    kappa2rel = (1-S1eff).^2 .* (1-S1eff.^((2+lambda_BC)./lambda_BC));
                elseif ismember(used_model,[7,8,9]) == true
                    % Douglas Jr. et al.
                    kappa1rel = (1 -  S1 ./(1-Sr2)).^2;
                    kappa2rel = ((S1-Sr1)./(1-Sr1)).^2;
                end

                omegac1 = eta1.*inv_T.*phi1 ./ (kappa.*rho1);
                omegac2 = eta2.*inv_T.*phi2 ./ (kappa.*rho2);
                kappa1  = kappa1rel .* kappa;
                kappa2  = kappa2rel .* kappa;

                % Implementing correction function
                if all(omegac1) == 0 && all(omegac2) == 0
                    C1 = ones(1,k);
                    C2 = ones(1,k);
                else
                    if eta1_in == 0 && eta2_in ~= 0
                        C1 = ones(1,k);
                        C2 = sqrt(1 + 0.5i .* omega ./ omegac2);
                    elseif eta1_in ~= 0 && eta2_in == 0
                        C1 = sqrt(1 + 0.5i .* omega ./ omegac1);
                        C2 = ones(1,k);
                    elseif eta1_in ~= 0 && eta2_in ~= 0
                        C1 = sqrt(1 + 0.5i .* omega ./ omegac1);
                        C2 = sqrt(1 + 0.5i .* omega ./ omegac2);
                    end
                end

                M             = (inv_N + phi1./K1 + phi2./K2 - ((S1.*S2.*phi) ./ (K1.*K2) + S1.*S2.*inv_N .* (S1 ./ K2 + S2 ./ K1)).*G).^(-1);
                M1            = M .* (1 - S1.*S2./K1 .* G);
                M2            = M .* (1 - S1.*S2./K2 .* G);
                Mtilde        = M .* (1 + inv_N.*S1.*S2./phi .* G);
                M1tilde       = M .* (1 - (inv_N.*S1.^2./phi + S1./K1) .* G);
                M2tilde       = M .* (1 - (inv_N.*S2.^2./phi + S2./K2) .* G);
                rho_star      = (1-phi).*rhos + (1-inv_T) .* (phi1.*rho1 + phi2.*rho2);
                gamma1        = phi1 .* inv_T;
                gamma2        = phi2 .* inv_T;
                Lambda1_star  = kappa1 ./ (C1.*phi1.*eta1);
                Lambda2_star  = kappa2 ./ (C2.*phi2.*eta2);
                Lambda11tilde = (rho_star .* Lambda1_star) ./ (phi1.*(inv_T.^2-inv_T));
                Lambda21tilde = (rho_star .* Lambda2_star) ./ (phi2.*(inv_T.^2-inv_T));
                Lambda12tilde = (inv_T ./ (rho1.*Lambda1_star) + 1./Lambda11tilde).^(-1);
                Lambda22tilde = (inv_T ./ (rho2.*Lambda2_star) + 1./Lambda21tilde).^(-1);
            end

            % Set attenuation, depending on porosity and viscosity
            attenuation = attenuation_default;
            if number_of_fluids == 0
                attenuation = false;
            elseif number_of_fluids == 1 && (phi_in == 0 || norm(eta1) == 0)
                attenuation = false;
            elseif number_of_fluids == 2
                if phi_in == 0 || (norm(eta1) == 0 && norm(eta2) == 0 )
                    attenuation = false;
                end
            end

            %% Wavespeed calculation
            % Allocate wavespeed  and inverse quality factor vectors
            v1 = zeros(k,1);
            v2 = zeros(k,1);
            v3 = zeros(k,1);

            Q1 = zeros(k,1);
            Q2 = zeros(k,1);
            Q3 = zeros(k,1);

            for  i = 1:k
                if number_of_fluids == 0
                    % Matrix elements of A
                    A11 = 0;
                    A12 = -(lambda_u(i)+2*my(i));
                    A21 = -1/rhos(i);
                    A22 = 0;

                    % Create matrices
                    A = [A11, A12; A21, A22];
                    E = zeros(2);

                elseif number_of_fluids == 1
                    if phi_in == 0
                        % Matrix elements of A
                        A11 = 0;
                        A12 = -(lambda_u(i)+2*my(i));
                        A13 = 0;
                        A14 = 0;

                        A21 = -1/rhos(i);
                        A22 = 0;
                        A23 = 0;
                        A24 = 0;

                        A31 = 0;
                        A32 = 0;
                        A33 = 0;
                        A34 = 0;

                        A41 = 1/rhos(i);
                        A42 = 0;
                        A43 = 0;
                        A44 = 0;

                        % Matrix E
                        E = zeros(4);
                    elseif phi_in ~= 0
                        % Matrix elements of A
                        A11 = 0;
                        A12 = -(lambda_u(i)+2*my(i)) + phi(i)*M(i)*b(i);
                        A13 = 0;
                        A14 = -phi(i)*M(i)*b(i);

                        A21 = -1/rho_star(i);
                        A22 = 0;
                        A23 = -gamma(i)/rho_star(i);
                        A24 = 0;

                        A31 = 0;
                        A32 = M(i)*(b(i)-phi(i));
                        A33 = 0;
                        A34 = phi(i)*M(i);

                        A41 = (inv_T(i) -1) / rho_star(i);
                        A42 = 0;
                        A43 = inv_T(i) / rho1(i) + (inv_T(i)*phi(i)*(inv_T(i)-1)) / rho_star(i);
                        A44 = 0;

                        % Matrix elements of E
                        E11 = 0;
                        E12 = 0;
                        E13 = 0;
                        E14 = 0;

                        E21 = 0;
                        E22 = -gamma(i)/(rho_star(i)*Lambda_star(i));
                        E23 = 0;
                        E24 = gamma(i)/(rho_star(i)*Lambda_star(i));

                        E31 = 0;
                        E32 = 0;
                        E33 = 0;
                        E34 = 0;

                        E41 = 0;
                        E42 = 1/Lambda2tilde(i);
                        E43 = 0;
                        E44 = -1/Lambda2tilde(i);
                    end

                    % Creating the matrices
                    A = [A11, A12, A13, A14; A21, A22, A23, A24; A31, A32, A33, A34; A41, A42, A43, A44];
                    if phi_in ~= 0
                        E = [E11, E12, E13, E14; E21, E22, E23, E24; E31, E32, E33, E34; E41, E42, E43, E44];
                    end

                elseif number_of_fluids == 2
                    if phi_in == 0
                        % Matrix elements of A
                        A11 = 0;
                        A12 = -(lambda_u(i)+2*my(i));
                        A13 = 0;
                        A14 = 0;
                        A15 = 0;
                        A16 = 0;

                        A21 = -1/rhos(i);
                        A22 = 0;
                        A23 = 0;
                        A24 = 0;
                        A25 = 0;
                        A26 = 0;

                        A31 = 0;
                        A32 = 0;
                        A33 = 0;
                        A34 = M(i)*S1(i)*S2(i)^2*G(i)/N(i);
                        A35 = 0;
                        A36 = -M(i)*S1(i)*S2(i)^2*G(i)/N(i);

                        A41 = 1/rhos(i);
                        A42 = 0;
                        A43 = 0;
                        A44 = 0;
                        A45 = 0;
                        A46 = 0;

                        A51 = 0;
                        A52 = 0;
                        A53 = 0;
                        A54 = -M(i)*S1(i)*S2(i)^2*G(i)/N(i);
                        A55 = 0;
                        A56 = M(i)*S1(i)*S2(i)^2*G(i)/N(i);

                        A61 = 1/rhos(i);
                        A62 = 0;
                        A63 = 0;
                        A64 = 0;
                        A65 = 0;
                        A66 = 0;

                        % Matrix E
                        E = zeros(6);
                    elseif phi_in ~= 0
                        % Matrix elements A
                        A11 = 0;
                        A12 = -(lambda_u(i)+2*my(i))+b(i)*(phi1(i)*M2(i)+phi2(i)*M1(i));
                        A13 = 0;
                        A14 = -phi1(i)*M2(i)*b(i);
                        A15 = 0;
                        A16 = -phi2(i)*M1(i)*b(i);

                        A21 = -1/rho_star(i);
                        A22 = 0;
                        A23 = -gamma1(i)/rho_star(i);
                        A24 = 0;
                        A25 = -gamma2(i)/rho_star(i);
                        A26 = 0;

                        A31 = 0;
                        A32 = M2(i)*b(i)-(phi1(i)*M2tilde(i)+phi2(i)*Mtilde(i));
                        A33 = 0;
                        A34 = phi1(i)*M2tilde(i);
                        A35 = 0;
                        A36 = phi2(i)*Mtilde(i);

                        A41 = (inv_T(i) -1) / rho_star(i);
                        A42 = 0;
                        A43 = inv_T(i)/rho1(i) + (inv_T(i)*phi1(i)*(inv_T(i)-1)) / rho_star(i);
                        A44 = 0;
                        A45 = (inv_T(i)*phi2(i)*(inv_T(i)-1)) / rho_star(i);
                        A46 = 0;

                        A51 = 0;
                        A52 = M1(i)*b(i)-(phi1(i)*Mtilde(i)+phi2(i)*M1tilde(i));
                        A53 = 0;
                        A54 = phi1(i)*Mtilde(i);
                        A55 = 0;
                        A56 = phi2(i)*M1tilde(i);

                        A61 = (inv_T(i) -1) / rho_star(i);
                        A62 = 0;
                        A63 = (inv_T(i)*phi1(i)*(inv_T(i)-1)) / rho_star(i);
                        A64 = 0;
                        A65 = inv_T(i)/rho2(i) + (inv_T(i)*phi2(i)*(inv_T(i)-1)) / rho_star(i);
                        A66 = 0;

                        % Matrix elements of E
                        E11 = 0;
                        E12 = 0;
                        E13 = 0;
                        E14 = 0;
                        E15 = 0;
                        E16 = 0;

                        E21 = 0;
                        E22 = -gamma1(i)/(rho_star(i)*Lambda1_star(i))-gamma2(i)/(rho_star(i)*Lambda2_star(i));;
                        E23 = 0;
                        E24 = gamma1(i)/(rho_star(i)*Lambda1_star(i));
                        E25 = 0;
                        E26 = gamma2(i)/(rho_star(i)*Lambda2_star(i));

                        E31 = 0;
                        E32 = 0;
                        E33 = 0;
                        E34 = 0;
                        E35 = 0;
                        E36 = 0;

                        E41 = 0;
                        E42 = 1/Lambda12tilde(i)+1/Lambda21tilde(i);
                        E43 = 0;
                        E44 = -1/Lambda12tilde(i);
                        E45 = 0;
                        E46 = -1/Lambda21tilde(i);

                        E51 = 0;
                        E52 = 0;
                        E53 = 0;
                        E54 = 0;
                        E55 = 0;
                        E56 = 0;

                        E61 = 0;
                        E62 = 1/Lambda11tilde(i) + 1/Lambda22tilde(i);
                        E63 = 0;
                        E64 = -1/Lambda11tilde(i);
                        E65 = 0;
                        E66 = -1/Lambda22tilde(i);
                    end

                    % Create the matrices
                    A = [A11, A12, A13, A14, A15, A16; A21, A22, A23, A24, A25, A26; A31, A32, A33, A34, A35, A36;
                         A41, A42, A43, A44, A45, A46; A51, A52, A53, A54, A55, A56; A61, A62, A63, A64, A65, A66];
                    if phi_in ~= 0
                        E = [E11, E12, E13, E14, E15, E16; E21, E22, E23, E24, E25, E26; E31, E32, E33, E34, E35, E36;
                             E41, E42, E43, E44, E45, E46; E51, E52, E53, E54, E55, E56; E61, E62, E63, E64, E65, E66];
                    end
                end

                %% Wavespeed calculation
                syms s % Wavenumber

                if number_of_fluids == 0
                    order_of_identity_matrix = 2;
                elseif number_of_fluids == 1
                    order_of_identity_matrix = 4;
                elseif number_of_fluids == 2
                    order_of_identity_matrix = 6;
                end

                E1       = 1i*E;                                                     % Matrix E gets its imaginary part
                equation = det( -omega(i)* eye(order_of_identity_matrix) + s*A +E1); % Solution of the similar eigenvalue problem
                solution = vpasolve(equation==0, s);                                 % Solves 'equation' with respect to the wavenumber s


                %% Calculation of v and Q
                v = omega(i)./real(solution);
                if number_of_fluids == 0 || phi_in == 0
                    P1    = find(v==max(v));
                    v1(i) = v(P1);

                    if attenuation == true
                        Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                    end
                elseif number_of_fluids == 1
                    P1    = find(v==max(v));
                    v1(i) = v(P1);
                    v(P1) = 0;
                    if length(solution) <= 2
                        v2(i) = 0;
                    else
                        P2    = find(v==max(v));
                        v2(i) = v(P2);
                    end

                    if attenuation == true
                        Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                        if length(solution) <= 2
                            Q2(i) = 0.
                        else
                            Q2(i) = 2*abs(imag(solution(P2))/real(solution(P2)));
                        end
                    end
                elseif number_of_fluids == 2
                    P1    = find(v==max(v));
                    v1(i) = v(P1);
                    v(P1) = 0;
                    if length(solution) <= 2
                        v2(i) = 0;
                    else
                        P2    = find(v==max(v));
                        v2(i) = v(P2);
                        v(P2) = 0;
                        if length(solution) > 4
                            P3    = find(v==max(v));
                            v3(i) = v(P3);
                        else
                            v3(i) = 0;
                        end
                    end

                    if attenuation == true
                        Q1(i) = 2*abs(imag(solution(P1))/real(solution(P1)));
                        if length(solution) <= 2
                            Q2(i) = 2*abs(imag(solution(P2))/real(solution(P2)));
                            if length(solution) > 4
                                Q2(i) = 2*abs(imag(solution(P3))/real(solution(P3)));
                            end
                        end
                    end
                end
            end

            if length(solution) == 2  % Definition of number_of_fluids depending on the porosity
                number_of_fluids = 0;
            elseif length(solution) == 4
                number_of_fluids = 1;
            elseif length(solution) == 6
                number_of_fluids = 2;
            end

            %% Calculation of the wavespeed for the most energetic angular frequency
            omega_f_c = linspace(min(omega),omega_p,2);

            if number_of_fluids == 0
                v1_f_c = interp1(omega,v1,omega_f_c);

                Q1_f_c = interp1(omega,Q1,omega_f_c);
            elseif number_of_fluids == 1
                v1_f_c = interp1(omega,v1,omega_f_c);
                v2_f_c = interp1(omega,v2,omega_f_c);

                Q1_f_c = interp1(omega,Q1,omega_f_c);
                Q2_f_c = interp1(omega,Q2,omega_f_c);
            elseif number_of_fluids == 2
                v1_f_c = interp1(omega,v1,omega_f_c);
                v2_f_c = interp1(omega,v2,omega_f_c);
                v3_f_c = interp1(omega,v3,omega_f_c);

                Q1_f_c = interp1(omega,Q1,omega_f_c);
                Q2_f_c = interp1(omega,Q2,omega_f_c);
                Q3_f_c = interp1(omega,Q3,omega_f_c);
            end

            %% GENERATING OUTPUT
            cd (workpath)
            if plot_diag
                run('plotting_output.m')
            end
            run('printing_output.m')

            %% Frequency stability
            %s            = min(v1) / (2.5*f_c);
            s            = min(v1./omega')*2*pi; % now: one element per wavelength
            n_cell_P1(u) = (max(mat_boundaries)-xmin) / s; % Frequency stability for P1

            if number_of_fluids >= 1
                %s            = min(v2) / (2.5*f_c);
                s            = min(v2./omega')*2*pi; % now: one element per wavelength
                n_cell_P2(u) = (max(mat_boundaries)-xmin) / s;
            end
            if number_of_fluids >= 2
                %s            = min(v3) / (2.5*f_c);
                s            = min(v3./omega')*2*pi; % now: one element per wavelength
                n_cell_P3(u) = (max(mat_boundaries)-xmin) / s;
            end
            n_cell(u) = max([n_cell_P1(u),n_cell_P2(u),n_cell_P3(u)]);
        end

        %% CALCULATE ELEMENT SIZES

        size_of_element = max(mat_boundaries)/max(n_cell); % size of one element
        fprintf('\n');
        fprintf('The minimal value for the number of elements to get frequency stability for each wave is %g.\n',ceil(max(n_cell)));
        if number_of_fluids >= 1
            fprintf('Frequency stability for P1 wave: %g \n',ceil(max(n_cell_P1)));
        end
        if number_of_fluids >= 2
            fprintf('Frequency stability for P2 wave: %g \n',ceil(max(n_cell_P2)));
        end
        fprintf('Size of one element to show each wave: %g \n',size_of_element);
        %fprintf('Frequency stability is calculated by the given most energetic frequency in "sources".\n');
        fprintf('\n');
    end
end

if logfile
diary off
end
cd (workpath)
cd ../tools
