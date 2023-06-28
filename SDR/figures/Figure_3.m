function Figure_3(gamma)

% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% This script generates Figure 3 for a given aspect ratio gamma
% To generate Figure 3 in the paper, use Figure_3(1.0)
%
% -----------------------------------------------------------------------------
% Authors: Matan Gavish and David Donoho <lastname>@stanford.edu, 2015
% 
% This program is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------


if nargin ==0
    gamma = 1.0;
end

loss_list = {{'F_1', 'F_2', 'F_3', 'F_4', 'F_5', 'F_6', 'F_7'},...
            {'O_1', 'O_2', 'O_3', 'O_4', 'O_5', 'O_6', 'O_7'}, ...
            {'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6', 'N_7'}, ... 
            {'Stein', 'Ent', 'Div', 'Fre', 'Aff'}} ;

names_list = {'Frobenius norm discrepancies', ...
              'Operator norm discrepancies', ...
              'Nuclear norm discrepancies', ...
              'Statistical discrepancies'};
jitter_factor = [1 2 1 1];

fs = 15;


h =  figure('Position',[100,100,1000,1300]);

for i_plot = 1:4
    axes_handle = subplot(2,2,i_plot);
    plot_shrinkers(gamma, loss_list{i_plot},axes_handle, ...
                        jitter_factor(i_plot) , 9)

    title(sprintf(' %s ($\\gamma=%0.2f$)', names_list{i_plot}, gamma), ...
                'FontSize',fs,'Interpreter','Latex');

end

figName = 'Figure_3';
print('-depsc2',[figName '.eps'])
saveas(h,[figName '.fig'], 'fig');

subs = ['f','o','n','s'];
for i_plot = 1:4
    h = figure;
    plot_shrinkers(gamma, loss_list{i_plot},gca, ...
                        jitter_factor(i_plot) , 9)

    title(sprintf(' %s ($\\gamma=%0.2f$)', names_list{i_plot}, gamma), ...
                'FontSize',fs,'Interpreter','Latex');
    
    figName = ['Figure_3' subs(i_plot)];
    print('-depsc2',[figName '.eps'])
    saveas(h,[figName '.fig'], 'fig');
end


