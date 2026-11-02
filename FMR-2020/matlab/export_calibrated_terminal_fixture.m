function export_calibrated_terminal_fixture(original_root,input_dir,output_dir)
% Reevaluate ORIGINAL HANKSOME solve_R.m at the already-calibrated Julia
% equilibrium prices, on the source 193/769 full asset grids. No fsolve.
if nargin~=3
    error('Use original_root, Julia input directory, output directory');
end
addpath(fullfile(original_root,'Steady_State_Ayagari'));
if ~exist('solve_R','file')
    error('Original Steady_State_Ayagari/solve_R.m is missing');
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
a=dlmread(fullfile(input_dir,'coarse_grid.csv'),',');
af=dlmread(fullfile(input_dir,'fine_grid.csv'),',');
ex=dlmread(fullfile(input_dir,'productivity.csv'),',');
P=dlmread(fullfile(input_dir,'transition.csv'),',');
v=dlmread(fullfile(input_dir,'prices.csv'),',');
params.beta=v(1); rate=v(2); terms=v(3); wage=v(4);
labor=v(5); profits=v(6);
params.gamma=v(7); params.theta=v(8); params.omega=v(9);
params.nA=length(a); params.nA_fine=length(af);
params.ns=length(ex);
params.Agrid=a; params.Agrid_fine=af;
params.bmin=a(33);
params.Index_b_min=33;
params.ex=ex; params.piex=P;
params.tol=1e-9; params.maxIter=5000;
assert(params.nA==193 && params.nA_fine==769 && params.ns==7);
initial_cpol=repmat(1+0.02*max(a,0),1,params.ns);
[~,cpol,~,cfine]=solve_R(rate,ex,labor,wage,terms,profits,initial_cpol,params);
dlmwrite(fullfile(output_dir,'source_coarse_policy.csv'),cpol,'precision','%.17g');
dlmwrite(fullfile(output_dir,'source_fine_policy.csv'),cfine,'precision','%.17g');
fprintf('AUTHOR CALIBRATED TERMINAL SOURCE POLICY rows=%d cols=%d\n',...
        size(cpol,1),size(cpol,2));
fprintf('Original source full-grid solve_R completed without 399-variable fsolve\n');
end
