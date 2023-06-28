
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
% 
%  compare analytic formula to numerical shrinker, for a specific value of
%  aspect ratio gamma, for all the analytic formulae provided in the paper.
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

addpath('../shrinkers')

available = {'F_1', 'F_2', 'F_3', 'F_4', 'F_6',...
             'O_1', 'O_2', 'O_6', ...
             'N_1', 'N_2', 'N_3', 'N_4', 'N_6', ...
             'Stein', 'Ent', 'Div', 'Fre', 'Aff'} ;

gamma = 1.0;

for loss = available 
    loss
    compare_numeric_analytic_shrinker(loss{1},gamma,1)
end
