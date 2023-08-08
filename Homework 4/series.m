%------------------------------------------%
%        *** SSSP - HOMEWORK #4 ***        %
%------------------------------------------%
%     3-port series scattering matrix      %
%------------------------------------------%
% Name: Iaccarino - Luan                   %
% Student ID: 10868500 - 10876787          %
%------------------------------------------%
function [seriesScatMat] = series(Z1, Z2, Z3)
seriesScatMat = zeros(3, 3);
ones = [1, 1, 1];
D = diag(ones);
vecZ = [Z1, Z2, Z3];
sumZ = sum(vecZ, "all");

seriesScatMat = D - (2/sumZ) * vecZ' * ones;
end