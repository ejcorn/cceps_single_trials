function pval = CIRC_CORRCL_P(alpha, x)

% wrapper of circ_corrcl that returns pval instead of rho as first output

[~,pval] = circ_corrcl_pairwise(alpha, x);
