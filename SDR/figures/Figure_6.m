function Figure_6(gamma)

% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% This script generates Figure 6 for a given aspect ratio gamma
% To generate Figure 6 in the supplemental article, use Figure_6(1.0)
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
   gamma=1;
end

lw3 = 1.3;
lfs = 13;
fs = 15;

table = tabulate_new(gamma);

table.PPI = (table.defaultLoss - table.bestLoss) ./ table.defaultLoss ;

h=figure;
colors    = colorDiscrepancies;
hold on
for i=1:length(table.loss_list),
 plot(table.lamList,100* table.PPI(:,i),colors{i},'LineWidth',lw3);
end
title(sprintf('PPI: Percentage Possible Improvement'),...
                  'FontSize',lfs,'Interpreter','Latex');

xlabel('$\lambda$','FontSize',lfs,'Interpreter','Latex');
ylabel('Improvement as a percentage of baseline',...
         'FontSize',lfs,'Interpreter','Latex');

xlim([3.8 10]);
ylim([0 100]);

figName = 'Figure_6a';
print('-depsc2',[figName '.eps'])
saveas(h,[figName '.fig'], 'fig');

h=figure;
colors    = colorDiscrepancies;
for i=1:length(table.loss_list),
 semilogx(table.lamList,100* table.PPI(:,i),colors{i},'LineWidth',lw3);
 if i==1
    hold on
   title(sprintf('PPI: Percentage Possible Improvement'),...
      'FontSize',lfs,'Interpreter','Latex');
   xlabel('$log(\lambda)$','FontSize',lfs,'Interpreter','Latex');
   ylabel('Improvement as a percentage of baseline',...
            'FontSize',lfs,'Interpreter','Latex');
   xlim([3.8 100]);
   ylim([0 100]);


 end
end

figName = 'Figure_6b';
print('-depsc2',[figName '.eps'])
saveas(h,[figName '.fig'], 'fig');

end

