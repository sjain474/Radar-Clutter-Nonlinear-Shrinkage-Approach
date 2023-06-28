
function table = tabulate_new(gamma)

% function table = tabulate_new(gamma)
%
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% Tabulate optimal shrinkers and their respective loss for all 26 loss functions 
% considered in the paper. This function is used to generate most of the figures
% in the paper and the supplemental article, but can be used regardless 
% to tabulate optimal shrinkers numerically.
% 
% IN: 
%   gamma: the desired aspect ratio
% 
% OUT:
%   table: a struct with fields 
%        table.loss_list : a list of the loss functions
%        table.loss_names  : a list of long names of the loss functions
%        table.lamList : a list of lambda values where shrinkers are calculated
%        table.bestLam  : the value of the optimal shrinker for loss i 
%                           at lambda value j
%        table.bestLoss : the loss of the optimal shrinkers as above
%        table.defaultLoss : the loss of hard thresholding for loss i
%                             at lambda value j, i.e no shrinkage
%        table.gamma : the value given as parameter
%
%
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



if gamma<=0 | gamma>1
    error('gamma must be between 0 and 1');
end 
if gamma*100 ~= round(gamma*100)
    error('please specify gamma to 2 decimal digits');
end

fname = sprintf('optimal_shrinkage_table_gamma=%0.2f.mat',gamma);



if exist(fname,'file') == 2
    fprintf('Found saved version in %s\n',fname);
    load(fname);
    return
else
    fprintf('no existing table found. making new table, this may take a while\n');
end

try
    matlabpool open
catch
end


% tabulate all optimal shrinkers and their respective losses for given gamma
loss_list = {'F_1', 'F_2', 'F_3', 'F_4', 'F_5', 'F_6', 'F_7',...
             'O_1', 'O_2', 'O_3', 'O_4', 'O_5', 'O_6', 'O_7', ...
             'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6', 'N_7', ... 
             'Stein', 'Ent', 'Div', 'Fre', 'Aff'} ;

loss_names = {'F+1', 'F+2', 'F+3', 'F+4', 'F+5', 'F+6', 'F+7',...
             'O+1', 'O+2', 'O+3', 'O+4', 'O+5', 'O+6', 'O+7', ...
             'N+1', 'N+2', 'N+3', 'N+4', 'N+5', 'N+6', 'N+7', ... 
             'Stein', 'Entropy', 'Divergence', 'Frechet', 'Affinity'} ;

%ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) + sqrt((lam+1-gamma).^2-4*lam))/(2.0));

lam_plus = (1+sqrt(gamma))^2;

lamList=unique([linspace(lam_plus,16,500) 16  20 30 40 50 60 ...
                    70 80 100 200 500 1000 2000 5000 10000]);
bestLam = zeros(length(lamList),length(loss_list));
bestLoss = zeros(length(lamList),length(loss_list));
defaultLoss = zeros(length(lamList),length(loss_list));

parfor i_loss = 1:length(loss_list)
    fprintf('loss=%s %d/%d\n',loss_names{i_loss},i_loss,length(loss_list));
    [bestLam(:,i_loss) bestLoss(:,i_loss) defaultLoss(:,i_loss)] = ...
                numerical_shrinkage(lamList, gamma,loss_list{i_loss});
end


table = struct();
table.loss_list = loss_list;
table.loss_names = loss_names;
table.lamList = lamList;
table.bestLam = bestLam;
table.bestLoss = bestLoss;
table.defaultLoss = defaultLoss;
table.gamma = gamma;

save(fname,'table');

end

