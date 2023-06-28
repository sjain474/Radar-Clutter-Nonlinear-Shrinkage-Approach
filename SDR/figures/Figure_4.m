function Figure_4(gamma)
    
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% This script generates Figure 4 for a given aspect ratio gamma
% To generate Figure 4 in the supplemental article, use Figure_4(1.0)
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


    if nargin == 0
        gamma = 1.0
    end

    if gamma<=0 | gamma>1
            error('gamma must be between 0 and 1');
    end 
    if gamma*100 ~= round(gamma*100)
        error('please specify gamma to 2 decimal digits');
    end
    
    lw3 = 1.2;
    lfs = 15;
    fs = 15;
    mw3 = 1.3;

    
    lam_plus = (1+sqrt(gamma))^2;
    ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) +  ...
                                        sqrt((lam+1-gamma).^2-4*lam))/(2.0));
 
    table = tabulate_new(gamma);
    
    h =  figure('Position',[100,100,1100,800]);
    colr = colorDiscrepancies;
    
    for jLoss=1:26
        plot(table.lamList,table.bestLam(:,jLoss),colr{jLoss},...
                    'LineWidth',lw3,'MarkerSize',mw3);
        if jLoss==1,
            xlabel('Empirical eigenvalue $\lambda$',...
                            'FontSize',lfs,'Interpreter','Latex');
            ylabel('Shrunken eigenvalue $\eta(\lambda)$',...
                        'FontSize',lfs,'Interpreter','Latex');
    
            title(sprintf(['Optimal Eigenvalue Shrinkers for Large Values ' ...
                        ' of $\\lambda$ ($\\gamma=%0.2f$)'],gamma),...
                               'FontSize',fs,'Interpreter','Latex');
            axis([0 100 0 100])
            hold on
        end
    end
    
    plot(table.lamList,table.lamList,'--','LineWidth',lw3);   % lambda
    plot(table.lamList,ell(table.lamList),'--','LineWidth',lw3);   % lambda
    h_legend = legend([table.loss_names(1:26) ...
                        {'$\eta(\lambda)=\lambda$' , '$\eta(\lambda)=\ell$'}]);
    set(h_legend,'Interpreter','Latex','FontSize',lfs,'Location','NorthWest');

    
    figName = 'Figure_4';
    print('-depsc2',[figName '.eps'])
    saveas(h,[figName '.fig'], 'fig');
    
end
