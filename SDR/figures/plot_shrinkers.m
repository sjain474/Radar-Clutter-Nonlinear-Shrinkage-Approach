function plot_shrinkers(gamma, loss_list,axes_handle,jitter_factor,max_lambda)

% function plot_shrinkers(gamma, loss_list,axes_handle,jitter_factor,max_lambda)
%
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% calculate optimal shrinkers (using numerical_shrinkage) and make a comparison
% plot. Used to generate Figures 2 and 3. You can use this function to 
% create comparison plots of your own.
% 
% IN:
% gamma: aspect ration to 2 decimal digits
% loss_list:  cell array of strings denoting loss functions whose corresponding 
%             optimal shrinkers are to be plotted. 
% axes_handle: handle to an active axes
% fitter_factor: how much to jitter in vertical axis. (optional, float).
%                0=no jitter 1=normal jitter
% max_lambda: plot from lambda_+ to this number on horizontal axis
%
% available loss functions: 
%             'F_1', 'F_2', 'F_3', 'F_4', 'F_5', 'F_6', 'F_7',
%             'O_1', 'O_2', 'O_3', 'O_4', 'O_5', 'O_6', 'O_7', 
%             'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6', 'N_7', 
%             'Stein', 'Ent', 'Div', 'Fre', 'Aff'
%
% plot is created without a title. add your title after running this function.
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

addpath('../shrinkers');
    
    if nargin < 4
        jitter_factor = 1;
    end
    if nargin <5
        max_lambda = 10;
    end
    
    if gamma<=0 | gamma>1
        error('gamma must be between 0 and 1');
    end 
    if gamma*100 ~= round(gamma*100)
        error('please specify gamma to 2 decimal digits');
    end
    
    table = tabulate_new(gamma);
    m = max(find(table.lamList <= max_lambda));
    
    
    % plot parameters
    fs = 15;
    lfs= 15;
    lfsl= 12;
    lw1 = 2;
    lw2 = 1.5;
    lw3 = 1.5;
    xlims = [3.5,max_lambda];
    ylims = [0.5,max_lambda];
    
    lam_plus = (1+sqrt(gamma))^2;
    ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) +  ...
                                        sqrt((lam+1-gamma).^2-4*lam))/(2.0));
    
    [found I] = ismember(loss_list,table.loss_list);
    if prod(found) ~= 1
        error('some loss names are not in table')
    end
    loss_list = table.loss_list(I);
    loss_names = table.loss_names(I);
    
    axes(axes_handle);
    hold on; 
    grid on;
    jitter = linspace(-0.10,0.10,length(I)) * jitter_factor;
    ColorSet=varycolor(length(I)+2);
    for i=1:length(I);
        plot(table.lamList(1:m), table.bestLam(1:m,I(i)) + jitter(i), ...
                         'Color', ColorSet(i,:),'LineWidth',lw1);
    
    
    end
    plot(table.lamList(1:m),table.lamList(1:m),'--',...
                'Color', ColorSet(length(I)+1,:),'LineWidth',lw2);   % lambda
    plot(table.lamList(1:m),ell(table.lamList(1:m)),'--',...
                'Color', ColorSet(length(I)+2,:),'LineWidth',lw2);   % ell
    xlim(xlims);
    ylim(ylims);
    
    xlabel('Empirical eigenvalue $\lambda$', ...
                                        'FontSize',lfs,'Interpreter','Latex');
    ylabel('Shrunken eigenvalue $\eta(\lambda)$', ...
                                        'FontSize',lfs,'Interpreter','Latex');
    h_legend=legend([loss_names ...
                    {'$\eta(\lambda)=\lambda$' , '$\eta(\lambda)=\ell(\lambda)$'}]);
    set(h_legend,'Interpreter','Latex','FontSize',lfsl,'Location','NorthWest');
    
end
