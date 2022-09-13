function lab = SCINOT_P(x)

% INPUTS:
% x: pp-value
% OUTPUTS:
% lab: number as text in scientific notation i.e. 2x10^-3

if x < 0.001
    lab = sprintf('%0.2f x 10^{%i}', 10^mod(log10(x),1),floor(log10(x)));
else
    lab = sprintf('%0.2g',x);
end