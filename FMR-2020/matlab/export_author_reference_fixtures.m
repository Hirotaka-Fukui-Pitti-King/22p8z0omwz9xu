function export_author_reference_fixtures(mat_dir, output_dir)
% Export the ORIGINAL author-solved MATLAB transition paths as full-precision
% CSV references consumed by scripts/evaluate_author_matlab_solution.jl.
% This function NEVER invokes fsolve. Run it where exp_flex.mat,
% exp_flex_lev16.mat, and exp_both_41_lev16.mat are available.
if nargin ~= 2
    error('Usage: export_author_reference_fixtures(mat_dir, output_dir)');
end
names = {'exp_flex.mat','exp_flex_lev16.mat','exp_both_41_lev16.mat'};
for k = 1:numel(names)
    if ~exist(fullfile(mat_dir,names{k}),'file')
        error('Missing author MATLAB solution: %s',fullfile(mat_dir,names{k}));
    end
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
base = load(fullfile(mat_dir,names{1}));
lev = load(fullfile(mat_dir,names{2}));
fix = load(fullfile(mat_dir,names{3}));
must_have(base,{'kl_path','G0_shock','params','kl_exp','resid','tot_exp'});
must_have(lev,{'kl_exp','resid','tot_exp','G0_shock'});
must_have(fix,{'k_l_t_fixed','resid','tot_fixed','G0_shock'});
TT = 200;
if numel(base.kl_path) ~= 2*TT-1 || size(base.G0_shock,1) ~= 769 || ...
        size(base.G0_shock,2) ~= 7
    error('Unexpected boom vector or shock distribution dimensions');
end
if max(abs(base.G0_shock(:)-lev.G0_shock(:))) > 1e-10 || ...
        max(abs(base.G0_shock(:)-fix.G0_shock(:))) > 1e-10
    error('Shock-date MATLAB distributions are inconsistent between runs');
end
write_array(output_dir,'original_boom_candidate',base.kl_path(:));
write_array(output_dir,'shock_pre_distribution',base.G0_shock);
p = base.params;
required = {'k_initial','q_initial','BG_shock','tot_initial', ...
            'k_final','beta','k_const_hat'};
for k = 1:numel(required)
    if ~isfield(p,required{k})
        error('Missing MATLAB params.%s',required{k});
    end
end
metadata = [p.k_initial;p.q_initial;p.BG_shock;p.tot_initial; ...
            p.k_final;p.beta;p.k_const_hat];
write_array(output_dir,'shock_inherited_metadata',metadata);
export_one(output_dir,'unleveraged',base.kl_exp,base.resid,base.tot_exp,2*TT-1,TT);
export_one(output_dir,'leveraged',lev.kl_exp,lev.resid,lev.tot_exp,2*TT-1,TT);
export_one(output_dir,'fixed',fix.k_l_t_fixed,fix.resid,fix.tot_fixed,3*TT-1,TT);
fprintf('AUTHOR MATLAB REFERENCE EXPORTED: %s (12 full-precision CSV files)\n',output_dir);
fprintf('No new GE solve performed; numerical residual checks use saved source outputs.\n');
end

function must_have(s,names)
for k=1:numel(names)
    if ~isfield(s,names{k})
        error('Saved MATLAB workspace lacks field: %s',names{k});
    end
end
end

function export_one(output_dir,label,x,residual,terms,dim,TT)
x=x(:); residual=residual(:); terms=terms(:);
if numel(x)~=dim || numel(residual)~=dim || numel(terms)~=TT
    error('Wrong source reference dimensions: %s',label);
end
if any(~isfinite(x)) || any(~isfinite(residual)) || ...
        any(~isfinite(terms)) || any(terms<=0)
    error('Nonfinite or nonpositive source reference: %s',label);
end
if max(abs(residual))>=1e-6
    error('Saved author vector is not a solved outer market system: %s',label);
end
write_array(output_dir,[label '_candidate'],x);
write_array(output_dir,[label '_source_residual'],residual);
write_array(output_dir,[label '_source_terms'],terms);
fprintf('%s: dimensions=%d, raw market max residual=%.17g, t1 terms=%.17g\n', ...
        label,dim,max(abs(residual)),terms(1));
end

function write_array(dirpath,name,values)
if any(~isfinite(values(:)))
    error('Nonfinite reference array: %s',name);
end
dlmwrite(fullfile(dirpath,[name '.csv']),values,'precision','%.17g');
end
