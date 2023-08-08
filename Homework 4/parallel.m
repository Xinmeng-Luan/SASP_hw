%------------------------------------------%
%        *** SSSP - HOMEWORK #4 ***        %
%------------------------------------------%
%    3-port parallel scattering matrix     %
%------------------------------------------%
% Name: Iaccarino - Luan                   %
% Student ID: 10868500 - 10876787          %
%------------------------------------------%
function [parallelScatMat] = parallel(Z1, Z2, Z3)
parallelScatMat = zeros(3, 3);
ones = [1, 1, 1];
D = diag(ones);
G1 = 1/Z1;
G2 = 1/Z2;
G3 = 1/Z3;
vecG = [G1, G2, G3];
sumG = sum(vecG, "all");

parallelScatMat = ((2/sumG) * vecG' * ones) - D;
parallelScatMat = parallelScatMat';
end