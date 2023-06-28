function Figure_5

% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% This script generates Figure 5 in the supplemental article
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




fname = 'OSAsympSlope.mat';

lw3 = 1.3;
lfs = 13;
fs = 15;

gammaList=linspace(.1,1,100);
bigLam = 100;
 
Discrep =   {'F_1', 'F_2', 'F_3', 'F_4', 'F_5', 'F_6', 'F_7', ...
             'O_1', 'O_2', 'O_3', 'O_4', 'O_5', 'O_6', 'O_7', ...
             'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6', 'N_7', ...
             'Stein', 'Ent', 'Div', 'Aff','Fre'} ;

if exist(fname,'file') == 2
    fprintf('Found saved version in %s\n',fname);
    load(fname);
    results = EmpiricalSlope;
else
    fprintf('no existing table found. making new table, this may take a while\n');
    
    results = zeros(length(gammaList),1+length(Discrep));
    
    for iGam= 1:length(gammaList),

        fprintf('%d/%d\n',iGam,length(gammaList));
        gamma = gammaList(iGam);
        lam_plus = (1+sqrt(gamma))^2;
        ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) + sqrt((lam+1-gamma).^2-4*lam))/(2.0));

        mu  = bigLam -1;
        results(iGam,1) = gamma;
        for iLoss=1:length(Discrep),
                eta = numerical_shrinkage(bigLam,gamma,Discrep{iLoss});
                asyShrink = eta/ell(bigLam);
                results(iGam,1+iLoss) = asyShrink;
        end
    end
    for iLoss=1:length(Discrep),
        fprintf('id=%i loss=%s asySlope(1)=%0.3f\n',iLoss,Discrep{iLoss},results(length(gammaList),1+iLoss));
    end
    
    EmpiricalSlope = results;
    save 'OSAsympSlope.mat' gammaList Discrep EmpiricalSlope
end    

h=figure;
colr = colorDiscrepancies;
hold on
for iLoss=1:length(Discrep),
    plot(results(:,1),results(:,1+iLoss),colr{iLoss},'LineWidth',lw3);
end
xlim([0.1 1]); ylim([0 1.01]);
xlabel('$\gamma$','FontSize',lfs,'Interpreter','Latex');
ylabel('asymptotic slope','FontSize',lfs,'Interpreter','Latex');
    
title(sprintf('Asymptotic Slopes of Optimal Shrinkers'),...
                               'FontSize',fs,'Interpreter','Latex');


figName = 'Figure_5';
print('-depsc2',[figName '.eps'])
saveas(h,[figName '.fig'], 'fig');

end
