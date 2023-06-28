function [bestLam bestLoss defaultLoss] = numerical_shrinkage(eigenvals,gamma,loss)

% function [bestLam bestLoss defaultLoss] = numerical_shrinkage(eigenvals,gamma,loss)
%
% See D. L. Donoho , M. Gavish and I. M. Johnstone,
% "OPTIMAL SHRINKAGE OF EIGENVALUES IN THE SPIKED COVARIANCE MODEL",
% https://arxiv.org/abs/1311.0851
%
% calculate optimal shrinker numerically. 
% 
% IN: 
%   eigenvals: a list of lambda values where the optimal shrinker should be
%     calculated.
%   gamma: aspect ratio
%   loss: the desired loss function
%
%    available loss functions: 
%             'F_1', 'F_2', 'F_3', 'F_4', 'F_5', 'F_6', 'F_7',
%             'O_1', 'O_2', 'O_3', 'O_4', 'O_5', 'O_6', 'O_7', 
%             'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6', 'N_7', 
%             'Stein', 'Ent', 'Div', 'Fre', 'Aff'
%
% OUT:
%  bestLam: the value of the optimal shrinker at the points given by eignevals
%  bestLoss: the corresponding optimal loss
%  defaultLoss: the corresponding loss of hard thresholding (no shrinkage) 
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


    bestLam = zeros(size(eigenvals));
    bestLoss = zeros(size(eigenvals));
    defaultLoss = zeros(size(eigenvals));
    for i=1:length(eigenvals)
        [bestLam(i) bestLoss(i) defaultLoss(i)] = ...
                numerical_shrinkage_impl(eigenvals(i),gamma,loss);
    end
end

function [bestLam bestLoss defaultLoss] = numerical_shrinkage_impl(lambda, gamma,loss)
    debug = 1;
    lam_plus = (1+sqrt(gamma))^2;

    if (lambda < lam_plus)
        error('lambda below bulk edge')
    end

    ell = @(lam) ((lam>=lam_plus).*((lam+1-gamma) + sqrt((lam+1-gamma).^2-4*lam))/(2.0));

    c = @(lam)((lam>=lam_plus) .* sqrt( (1-gamma./((ell(lam)-1).^2) ) ./ (1+gamma./(ell(lam)-1) )   ));
 
        
    half = @(A,B)( inv(sqrt(A)) * B * inv(sqrt(A)) );
    Delta_1 = @(A,B)(A-B);
    Delta_2 = @(A,B)(inv(A)-inv(B));
    Delta_3 = @(A,B)(inv(A)*B-eye(size(A)));
    Delta_4 = @(A,B)(inv(B)*A-eye(size(A)));
    Delta_5 = @(A,B)(inv(A)*B + inv(B)*A-2*eye(size(A)));
    Delta_6 = @(A,B)(half(A,B) - eye(size(A)));
    Delta_7 = @(A,B)(logm(half(A,B)));
    
    impl_F_1 = @(A,B)(norm( Delta_1(A,B) ,'fro'));
    impl_F_2 = @(A,B)(norm( Delta_2(A,B) ,'fro'));
    impl_F_3 = @(A,B)(norm( Delta_3(A,B) ,'fro'));
    impl_F_4 = @(A,B)(norm( Delta_4(A,B) ,'fro'));
    impl_F_5 = @(A,B)(norm( Delta_5(A,B) ,'fro'));
    impl_F_6 = @(A,B)(norm( Delta_6(A,B) ,'fro'));
    impl_F_7 = @(A,B)(norm( Delta_7(A,B) ,'fro'));
    
    
    impl_O_1 = @(A,B)(norm( Delta_1(A,B) ));
    impl_O_2 = @(A,B)(norm( Delta_2(A,B) ));
    impl_O_3 = @(A,B)(norm( Delta_3(A,B) ));
    impl_O_4 = @(A,B)(norm( Delta_4(A,B) ));
    impl_O_5 = @(A,B)(norm( Delta_5(A,B) ));
    impl_O_6 = @(A,B)(norm( Delta_6(A,B) ));
    impl_O_7 = @(A,B)(norm( Delta_7(A,B) ));
 
    impl_N_1 = @(A,B)(norm( svd(Delta_1(A,B)),1 ));
    impl_N_2 = @(A,B)(norm( svd(Delta_2(A,B)),1 ));
    impl_N_3 = @(A,B)(norm( svd(Delta_3(A,B)),1 ));
    impl_N_4 = @(A,B)(norm( svd(Delta_4(A,B)),1 ));
    impl_N_5 = @(A,B)(norm( svd(Delta_5(A,B)),1 ));
    impl_N_6 = @(A,B)(norm( svd(Delta_6(A,B)),1 ));
    impl_N_7 = @(A,B)(norm( svd(Delta_7(A,B)),1 ));

    impl_Stein =  @(A,B)(SteinLoss(A,B));
    impl_Ent = @(A,B)(SteinLoss(B,A));
    impl_Div = @(A,B)( trace(inv(A) * B  + inv(B) *A - 2*eye(size(A)))/2.0 );
    impl_Fre = @(A,B)(trace(A + B - 2.0* sqrtm(A)*sqrtm(B)));
    impl_Aff = @(A,B)(0.5*log(0.5*det(A+B)/(sqrt(det(A))*sqrt(det(B)))));
    
    A = [ell(lambda) 0 ; 0 1];

    eval(['LossFunc = @(A,B)(impl_' loss '(A,B));']);

    %result = bestLam(lambda,ell(lambda),A,c(lambda),LossFunc,debug);
    [bestLam bestLoss] = bestLam_impl(A,c(lambda),LossFunc);
    defaultLoss = defaultLoss_impl(lambda,A,c(lambda),LossFunc);
    
end

function J = SteinLoss(A,B)
Delta = sqrtm(A)^(-1) * B * sqrtm(A)^(-1);
J = (trace(Delta) - log(det(Delta)) - trace(eye(size(A))))./2;
end

function [optLam optVal] = bestLam_impl(A,c,J)
    lobnd=1;
    upbnd=A(1,1)+2;
    s = -sqrt(1 - c^2);
    u = [c;s];
    optVal = Inf;
    optLam = 0;
   for iter=1:6,
    lamList = linspace(lobnd,upbnd,100);
    inx     = 0;
    for iLam=1:length(lamList),
        lam = lamList(iLam);
        eta = lam-1;
        B = eye(size(A)) + eta*(u*u');
        val = J(A,B);
        if val < optVal,
            optVal = val;
            optLam = lam;
            inx = iLam;
        end
    end
    if inx>0,
        loinx=max(1,inx-1);
        hiinx=min(inx+1,length(lamList));
        lobnd=lamList(loinx);
        upbnd=lamList(hiinx);
    else
        break;
    end
   end

end

function risk = defaultLoss_impl(lam,A,c,J)
    
    s = -sqrt(1 - c^2);
    u = [c;s];

    A1   = eye(size(A)) + (lam-1) .* u*(u');
    risk = J(A,A1);
end
