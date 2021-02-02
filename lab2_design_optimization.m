%% Solar Cell Front Contact Grid Design
%
%  Done for ELEC4703, Design Lab 2
%  Instructor: Ravi Prakash
%
%  Date: 2nd Feb 2021
%  Author: Ragini Bakshi

clear all;
close all;

%% Defining Parameters
Rss = 60;               % given sheet resistance in Ohm/square
lmda = 4e-6;            % design dimension defined by fabrication in m
cellW = 10.16e-3/lmda;  % cell active area width in m / lambda
ag_t = 1.2e-4;          % silver thickness in cm
ag_r = 1.6e-6;          % silver resistivity in ohm*cm
rho_s = ag_r/ag_t;      % metal sheet resistance (for busbar and fingers)

J_mp = 0.03519*10^4;     % obtained from PC1D simulator in A/m^2
V_mp = 0.53;             % obtained from PC1D simulator in V

% Determine:
%S = 0;                  % spacing b/w fingers (center to center)
%optimized_B = 0;                  % finger length: set to S/2 less than full width of cell
%optimized_A = 0;                  % busbar length = width of half cell
%n_f --> number of fingers
%Wf = 0;                 % finger width
%Wb = 0;                 % half width of busbar

%% Simulation Settings
% Horizontal Layout : 0
% Vertical Layout: 1
busbar = 0;

% Finger Geometry: rectangular = 3 vs triangular = 4
m = 3;
% Busbar Geometry: rectangular = 3 vs triangular = 4 
n = 3;

%%
min_loss = 1000;
for Wf = lmda : lmda : 10*lmda
    for Wb = lmda : lmda : cellW/4*lmda
        for s = lmda : lmda : cellW/2*lmda
            if busbar == 1              % vertical busbar
                A = cellW*lmda;
                B = cellW*lmda/2 - s/2;
               
            else                        % horizontal busbar
                A = cellW*lmda/2;       % half the cell width since busbar is horizontal
                B = cellW*lmda - s/2;   % rule of thumb
            end
            
            p_re = (Rss*s^2*J_mp)/(12*V_mp);            % power loss due to current flow thru emitter layer resistance to grid
            p_sf = Wf/s;                                % power loss due to shadowing by grid fingers
            p_rf = (B^2*rho_s*J_mp*s)/(m*V_mp*Wf);      % power loss due to current flow thru grid finger resistance
            p_sb = Wb/B;                                % power loss due to shadowing by busbar
            p_rb = (A^2*B*rho_s*J_mp)/(n*V_mp*Wb);      % power loss due to current flow thur busbar resistance
            
            total_losses = p_re + p_sf + p_rf + p_sb + p_rb;
            
            if(total_losses < min_loss)
                s_new = s;
                B_new = B;
                A_new = A;
                Wf_new = Wf;
                Wb_new = Wb;
                min_loss = total_losses;
                p_re_new = p_re;
                p_sf_new = p_sf;
                p_rf_new = p_rf;
                p_sb_new = p_sb;
                p_rb_new = p_rb;
            end
        end
    end
end

%% Optimized Parameters in terms of Lambda for L-Edit:
optimized_B = B_new / lmda;
optimized_A = A_new / lmda;
optimized_S = s_new / lmda;
optimized_Wf = Wf_new / lmda;
optimized_Wb = Wb_new / lmda;

n_f = cellW/2/optimized_S;       % number of fingers
table_1 = table(Rss, n_f, optimized_B, optimized_A, optimized_S, optimized_Wf, optimized_Wb, total_losses)
