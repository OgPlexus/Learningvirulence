function tex = makeAltMethodsLatexTable(A, outFile, caseLabel, tabLabel)
% makeAltMethodsLatexTable  Write a LaTeX table comparing GeoC with transfer entropy, CCM and linear Granger.
%
% CALLS     None
% CALLED BY runCase1a
%
% INPUT
%   A          struct with nSurrogates, Lmin, Lmax, geoc.b2m, geoc.m2b (scalars), te.b2m and te.m2b
%              (transferEntropyKSG outputs), ccm.b2m and ccm.m2b (ccmSkill outputs), gr.b2m and gr.m2b
%              (grangerLinearF outputs). b2m is beta to mu and m2b is mu to beta.
%   outFile    path of the .tex file
%   caseLabel  char used in the caption, for example 'Case~1(a)'
%   tabLabel   char LaTeX label of the table
% OUTPUT
%   tex        table text, also written to outFile
pf = @(p) ternary(p < 0.001, '$p<0.001$', sprintf('$p=%.3f$', p));
L = {};
L{end+1} = '\begin{table}[htbp]';
L{end+1} = '\centering';
L{end+1} = '\small';
L{end+1} = sprintf(['\\caption{Comparison of directional measures for %s, in which virulence influences transmission ' ...
  'by construction and transmission does not enter the virulence update. Each entry is computed from the same ' ...
  'simulated record. GeoC is the raw directional estimate. Transfer entropy (TE) is a $k$-nearest-neighbor estimate ' ...
  'in nats with a source-permutation $p$ value (%d surrogates). For convergent cross mapping (CCM), $\\rho$ is the ' ...
  'cross-map skill at the largest library size $L=%d$ and $\\Delta\\rho$ is its increase from $L=%d$. ' ...
  'The linear Granger test is an $F$ test with one lag.}'], caseLabel, A.nSurrogates, A.Lmax, A.Lmin);
L{end+1} = sprintf('\\label{%s}', tabLabel);
L{end+1} = '\begin{tabular}{llcc}';
L{end+1} = '\toprule';
L{end+1} = 'Method & Statistic & $\beta\to\mu$ (no pathway) & $\mu\to\beta$ (pathway) \\';
L{end+1} = '\midrule';
L{end+1} = sprintf('GeoC & $\\mathrm{GeoC}_{Y\\to X}$ & $%.4f$ & $%.4f$ \\\\', A.geoc.b2m, A.geoc.m2b);
L{end+1} = sprintf('Transfer entropy & TE (nats), %s & $%.3f$ (%s) & $%.3f$ (%s) \\\\', 'permutation $p$', ...
    A.te.b2m.TE, pf(A.te.b2m.p), A.te.m2b.TE, pf(A.te.m2b.p));
L{end+1} = sprintf('CCM & $\\rho$ ($\\Delta\\rho$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) \\\\', ...
    A.ccm.b2m.rho(end), A.ccm.b2m.deltaRho, A.ccm.m2b.rho(end), A.ccm.m2b.deltaRho);
L{end+1} = sprintf('Linear Granger & $F$ (%s) & $%.1f$ (%s) & $%.1f$ (%s) \\\\', 'p value', ...
    A.gr.b2m.F, pf(A.gr.b2m.p), A.gr.m2b.F, pf(A.gr.m2b.p));
L{end+1} = '\bottomrule';
L{end+1} = '\end{tabular}';
L{end+1} = '\end{table}';
tex = strjoin(L, newline);
fid = fopen(outFile, 'w'); fprintf(fid, '%s\n', tex); fclose(fid);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
