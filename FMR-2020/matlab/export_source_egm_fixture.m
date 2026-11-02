function export_source_egm_fixture(original_root, output_dir)
% Original-MATLAB-source one-date fixture, no 399-variable root solve.
% Run from MATLAB, with original_root pointing at kurtmitman/HANKSOME.
% The source's Rouwenhorst implementation is used, NOT a second port.
if nargin < 1
    error('Pass the path to the original kurtmitman/HANKSOME repository');
end
if nargin < 2
    output_dir = fullfile(fileparts(mfilename('fullpath')), '..', ...
        'results', 'matlab_egm_fixture');
end
addpath(fullfile(original_root,'Exogenous_files'));
addpath(fullfile(original_root,'Steady_State_Ayagari'));
if ~exist('rouwenhorst','file')
    error('Original Exogenous_files/rouwenhorst.m not found');
end
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

alpha=0.33;
gamma=1;
rho=0.97;
sigma_y=0.84*sqrt(1-rho^2);
[log_ex,P]=rouwenhorst(7,0,rho,sigma_y);
ex=exp(log_ex);
stationary=P^1000000;
ex=ex/(stationary(1,:)*ex);

gdp=12.65/4;
nfa=-6.5*0.69;
capital_pc=26.2-nfa;
ky=capital_pc/gdp;
Kref=ky^(1/(1-alpha));
bmin=-1e-7*Kref;
amax=999*Kref;
z_coarse=linspace(-0.2,1,193)';
z_fine=linspace(-0.2,1,769)';
z_coarse(33)=0;
z_fine(129)=0;
Agrid=bmin+(amax-bmin)*z_coarse.^3;
Agrid_fine=bmin+(amax-bmin)*z_fine.^3;

% One-step test inputs also documented in the Julia comparison script.
beta=.98;
r_now=.01;
r_next=.01;
terms=1;
wage=1.8;
labor=1;
profits=.25;
omega=.4;
theta=1;
consumption_floor=1e-5;
% Positive and increasing with assets: a_today spans the fine/coarse
% source grid so interp1(...,'pchip') is tested WITHIN its domain.
% A constant future policy caused an unrelated source upper-tail
% out-of-domain NaN, not a valid PCHIP equivalence experiment.
future=repmat(1+0.02*max(Agrid,0),1,7);
m=1+omega/(1-omega)*terms^(theta-1);

coarse_consumption=zeros(193,7);
fine_consumption=zeros(769,7);
fine_savings=zeros(769,7);
fine_budget=zeros(769,7);
for j=1:7
    income=wage*labor*ex(j)+profits;
    emu=future.^(-gamma)*P(j,:)';
    c_s=(beta*(1+r_next)*emu).^(-1/gamma);
    a_today=(m*c_s+Agrid-income)/(1+r_now);
    c_constrained=max(((1+r_now)*Agrid+income-bmin)/m,consumption_floor);
    threshold=a_today(33);
    for i=1:193
        if Agrid(i)<=threshold
            coarse_consumption(i,j)=c_constrained(i);
        else
            coarse_consumption(i,j)=max( ...
                interp1(a_today,c_s,Agrid(i),'pchip'),consumption_floor);
        end
    end
    for i=1:769
        if Agrid_fine(i)<=threshold
            fine_consumption(i,j)=interp1(Agrid,c_constrained, ...
                                         Agrid_fine(i),'pchip');
            fine_savings(i,j)=bmin;
        else
            fine_consumption(i,j)=interp1(a_today,c_s, ...
                                         Agrid_fine(i),'pchip');
            fine_savings(i,j)=(1+r_now)*Agrid_fine(i)+income- ...
                              m*fine_consumption(i,j);
        end
        fine_budget(i,j)=(1+r_now)*Agrid_fine(i)+income- ...
                           m*fine_consumption(i,j)-fine_savings(i,j);
    end
end
% No estimate of an equilibrium, only matched one-step author formulas.
dlmwrite(fullfile(output_dir,'coarse_grid.csv'),Agrid,'precision','%.17g');
dlmwrite(fullfile(output_dir,'fine_grid.csv'),Agrid_fine,'precision','%.17g');
dlmwrite(fullfile(output_dir,'productivity.csv'),ex,'precision','%.17g');
dlmwrite(fullfile(output_dir,'transition.csv'),P,'precision','%.17g');
dlmwrite(fullfile(output_dir,'coarse_consumption.csv'),coarse_consumption,'precision','%.17g');
dlmwrite(fullfile(output_dir,'fine_consumption.csv'),fine_consumption,'precision','%.17g');
dlmwrite(fullfile(output_dir,'fine_savings.csv'),fine_savings,'precision','%.17g');
dlmwrite(fullfile(output_dir,'fine_budget.csv'),fine_budget,'precision','%.17g');

% Independently test the full-grid stationary policy produced by the
% original solve_R.m. Same positive increasing initial c0 as Julia.
params.nA=193; params.nA_fine=769; params.ns=7;
params.Agrid=Agrid; params.Agrid_fine=Agrid_fine;
params.Index_b_min=33;
params.bmin=bmin;
params.ex=ex; params.piex=P;
params.gamma=gamma; params.theta=theta; params.omega=omega;
params.beta=beta; params.tol=1e-9; params.maxIter=5000;
initial_cpol=future;
[~,cpol_source,~,cfine_source]=solve_R( ...
    r_now,ex,labor,wage,terms,profits,initial_cpol,params);
dlmwrite(fullfile(output_dir,'stationary_coarse_consumption.csv'), ...
         cpol_source,'precision','%.17g');
dlmwrite(fullfile(output_dir,'stationary_fine_consumption.csv'), ...
         cfine_source,'precision','%.17g');
fprintf('MATLAB source solve_R stationary fixture complete\n');
fprintf('MATLAB EGM FIXTURE SAVED: %s\n',output_dir);
fprintf('GRID: %d coarse, %d fine, constraint indices 33 and 129\n', ...
        length(Agrid),length(Agrid_fine));
end
