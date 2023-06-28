function result = compare_numeric_analytic_shrinker_(loss,gamma,show)

% function result = compare_numeric_analytic_shrinker_(loss,gamma,show)
%
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% compare the optimal shrinker from the analytic formula to a numeric
% computation of the optimal shrinker
%
% IN:
%  loss: a string indicating loss function. one of -
%
%             F_1: Frobenius norm on A-B
%             F_2: Frobenius norm on precision (A^-1 - B^-1)
%             F_3: Frobenius norm on A^-1 * B - I
%             F_4: Frobenius norm on B^-1 * A - I
%             F_6 Frobenius norm on A^1/2 * B * A^1/2 - I
%             N_1: Nuclear norm on A-B
%             N_2: Nuclear norm on precision (A^-1 - B^-1)
%             N_3: Nuclear norm on A^-1 * B - I
%             N_4: Nuclear norm on B^-1 * A - I
%             N_6 Nuclear norm on A^1/2 * B * A^1/2 - I
%             O_1: Operator norm on A-B
%             O_2: Operator norm on precision (A^-1 - B^-1)
%             O_6 Operator norm on A^1/2 * B * A^1/2 - I
%             Stein: Stein's loss
%             Ent: Entropy loss
%             Div: Divergence loss
%             Fre: Frechet loss
%             Aff: Affine loss
%
%  lam_vals: a vector of lambda values to compare. 
%  gamma: aspect ratio
%  show: 0-no figure, 1-comparison figure
% 
% OUT:
%   result: the l2 difference between the numerical and analytic optimal
%   shrinkers, per entry.
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

    if nargin < 3
        show = 0;
    end

    lam_plus = (1+sqrt(gamma))^2;
    num_vals = 200;
    lam_vals = linspace(lam_plus,4*lam_plus,num_vals);
    result = compare_numeric_analytic_shrinker_impl(loss,lam_vals,gamma,show);
end


function result= compare_numeric_analytic_shrinker_impl(loss,lam_vals,gamma,show)


available = {'F_1', 'F_2', 'F_3', 'F_4', 'F_6', 'O_1', 'O_2', 'O_6', 'N_1', 'N_2', ... 
'N_3', 'N_4', 'N_6', 'Stein', 'Ent', 'Div', 'Fre', 'Aff'} ;

if ismember(loss,available)

    analytic_eta = optimal_shrinkage(lam_vals,gamma,loss,1);
else
    analytic_eta = zeros(size(lam_vals));
end

    
    numerical_eta = numerical_shrinkage(lam_vals, gamma, loss);
    
    if(show) 
        figure;
        plot(lam_vals,analytic_eta,'ob',lam_vals,numerical_eta,'+r');
        title(loss)
    end
    assert(length(analytic_eta)==length(numerical_eta));
    result = norm(analytic_eta - numerical_eta) / length(numerical_eta); 
end
