function export_original_399_residual(original_root,input_dir,terminal_dir,output_dir,mode)
% Evaluate the UNMODIFIED original author's 399-equation residual function
% at precisely the K/L path and shock-date distribution exported by Julia.
% Do not call outer fsolve. Unscaled original source units are preserved.
if nargin<4 || nargin>5
    error('Expected original_root, input_dir, terminal_dir, output_dir, [mode]');
end
if nargin==4, mode='literal'; end
if ~strcmp(mode,'literal') && ~strcmp(mode,'tail_extrapolation')
    error('Unsupported original MATLAB residual evaluation mode');
end
addpath(fullfile(original_root,'Transition'));
addpath(fullfile(original_root,'Exogenous_files'));
% Octave has fsolve but lacks MATLAB's optimoptions. This compatibility
% adapter only translates the scalar fsolve options; the ORIGINAL author
% economics and residual function files are unchanged.
if exist('OCTAVE_VERSION','builtin')
    addpath(fullfile(fileparts(mfilename('fullpath')),'octave_compat'));
end
if strcmp(mode,'tail_extrapolation')
    % Source-compatibility TEST variant, not the original: avoids only
    % out-of-support PCHIP NaN at fine-grid nodes. It MUST NOT be
    % described as literal author execution, and equality is not assumed.
    addpath(fullfile(fileparts(mfilename('fullpath')),
                     'octave_compat','source_support'));
    fprintf('REFERENCE VARIANT ONLY: fine-grid PCHIP tail extrapolation enabled\n');
else
    fprintf('LITERAL ORIGINAL backsolve_egm.m with no source patches\n');
end
if ~exist('solve_transition_bis','file') || ~exist('basefun','file')
    error('Original source Transition/solve_transition_bis.m or basefun.m missing');
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
x=dlmread(fullfile(input_dir,'candidate_399.csv'),',');
B=dlmread(fullfile(input_dir,'B_stop.csv'),',');
z=dlmread(fullfile(input_dir,'z_path.csv'),',');
a=dlmread(fullfile(input_dir,'coarse_grid.csv'),',');
af=dlmread(fullfile(input_dir,'fine_grid.csv'),',');
G0=dlmread(fullfile(input_dir,'mu_pre.csv'),',');
ex=dlmread(fullfile(input_dir,'productivity.csv'),',');
P=dlmread(fullfile(input_dir,'transition.csv'),',');
v=dlmread(fullfile(input_dir,'source_parameters.csv'),',');
meta=dlmread(fullfile(input_dir,'metadata.csv'),',');
cpol=dlmread(fullfile(terminal_dir,'source_coarse_policy.csv'),',');
cfine=dlmread(fullfile(terminal_dir,'source_fine_policy.csv'),',');
TT=length(B);
assert(TT==200 && numel(x)==2*TT-1);
assert(numel(z)==TT && size(G0,1)==769 && size(G0,2)==7);
assert(numel(a)==193 && numel(af)==769 && numel(ex)==7);
assert(size(P,1)==7 && size(P,2)==7);
assert(size(cpol,1)==193 && size(cpol,2)==7);
assert(size(cfine,1)==769 && size(cfine,2)==7);
assert(numel(v)==15 && size(meta,1)==2 && size(meta,2)==11);
params.TT=TT;
params.beta=v(1);
params.omega=1-v(2);
params.theta=v(3);
params.gamma=v(4);
params.alpha=v(5);
params.delta=v(6);
params.epsilon_f=v(7);
params.epsilon_w=v(8);
params.phi_ac=v(9);
params.theta_star=v(10);
params.d_s=v(11);
params.chi_dis=v(12);
params.k_final=v(13);
params.linsp_min=v(14);
params.curv_fine=v(15);
params.Agrid=a(:);
params.Agrid_fine=af(:);
params.Amax=af(end);
params.bmin=a(33);
params.Index_b_min=33;
params.nA=numel(a);
params.nA_fine=numel(af);
params.ns=numel(ex);
params.nx=numel(ex);
params.ex=ex(:);
params.piex=P;
params.ex_mean=1;
params.k_initial=meta(1,1);
params.q_initial=meta(1,2);
params.BG_shock=meta(1,3);
params.tot_initial=meta(1,4);
params.k_const_hat=0;
params.i_ss=params.delta*params.k_final;
params.b_agg_ss=B(end);
params.fixed=false;
params.foreign=true;
params.theta_nr=100;
params.tol_errore_t=1e-6;
% This is the original MATLAB implementation without reinterpretation:
% original scalar fsolve(x0=.8), full 193-node EGM, full 769-node
% distribution, original c>=1e-5 floor and author residuals.
fprintf('ORIGINAL 399 SOURCE DIRECT EVALUATION begin, T=%d\n',TT);
[resid,rk,rate,aprime,assets,Y,cbudget,netc,tot,ct, ...
 lab,indis,wage,dcoarse,dfine,cfine_tran,profits,q,dist] = ...
 solve_transition_bis(x(:),z(:),B(:),cpol,cfine,G0,params);
resid=resid(:);
dlmwrite(fullfile(output_dir,'matlab_raw_399.csv'),resid, ...
         'precision','%.17g');
dlmwrite(fullfile(output_dir,'matlab_terms.csv'),tot(:), ...
         'precision','%.17g');
dlmwrite(fullfile(output_dir,'matlab_aggregate_C.csv'),ct(:), ...
         'precision','%.17g');
dlmwrite(fullfile(output_dir,'matlab_aggregate_A.csv'),aprime(:), ...
         'precision','%.17g');
dlmwrite(fullfile(output_dir,'matlab_rate.csv'),rate(:), ...
         'precision','%.17g');
dlmwrite(fullfile(output_dir,'matlab_wage.csv'),wage(:), ...
         'precision','%.17g');
if any(~isfinite(resid))
    fprintf('ORIGINAL SOURCE 399 RESIDUAL UNDEFINED: %d of %d nonfinite; do not claim parity or equilibrium\n', ...
            sum(~isfinite(resid)),length(resid));
else
    fprintf('ORIGINAL SOURCE 399 RAW RESIDUAL FINITE, max=%g norm=%g\n', ...
            max(abs(resid)),norm(resid));
end
end
