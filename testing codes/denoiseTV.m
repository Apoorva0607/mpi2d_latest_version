function [x,J] = denoiseTV(y,lambda,Nit)
% [x,J] = denoiseTVclip(y,lambda,a,Nit)
% Total variation filtering (denoising)
% INPUT
%   y - noisy signal (row vector)
%   lambda - regularization parameter
%   Nit - Number of iterations
% OUTPUT
%   x - result of denoising
%   J - objective function

J = zeros(1,Nit);       % Objective function
N = length(y);
z = zeros(1,N-1);       % Initialize z
alpha = 4;
for k = 1:Nit
    Dtz = [-z(1) -diff(z) z(end)];    % D'*z
    x = y - (lambda/2) * Dtz;
    J(k) = sum(abs(x-y).^2) + lambda * sum(abs(diff(x)));
    z = z + 2/(lambda*alpha) * diff(x);
    
end