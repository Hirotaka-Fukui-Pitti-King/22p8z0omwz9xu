function export_author_fixed_audit(mat_dir, output_dir)
% Extract exact ORIGINAL fixed-rate time paths and selected policy arrays.
% Only reads author-saved .mat; never re-solves any equilibrium.
if nargin~=2
    error('Expected mat_dir and output_dir');
end
filename=fullfile(mat_dir,'exp_both_41_lev16.mat');
if ~exist(filename,'file'), error('Original MAT missing: %s',filename); end
s=load(filename);
needed={'k_l_t_fixed','resid','r_t_fixed','w_nr_guess_fixed', ...
        'Profits_fixed','q_agg_fixed','In_dis_fixed', ...
        'c_tran_fine_fixed','dec_tran_fine_fixed','tot_fixed', ...
        'new_shock_b_agg','cpol','c_fine'};
for i=1:numel(needed)
    if ~isfield(s,needed{i})
        error('Fixed MATLAB saved workspace missing variable %s',needed{i});
    end
end
T=200;
assert(numel(s.k_l_t_fixed)==3*T-1 && numel(s.resid)==3*T-1);
assert(size(s.In_dis_fixed,1)==769 && size(s.In_dis_fixed,2)==7);
assert(size(s.c_tran_fine_fixed,1)==769 && size(s.c_tran_fine_fixed,2)==7);
assert(size(s.dec_tran_fine_fixed,1)==769 && size(s.dec_tran_fine_fixed,2)==7);
assert(size(s.c_tran_fine_fixed,3)==T);
% Original backsolve_egm fills savings decisions only for t=1,...,T-1.
% Its terminal consumption is supplied externally; no terminal savings policy
% is constructed. The saved MATLAB savings array therefore has T-1 slices.
assert(size(s.dec_tran_fine_fixed,3)==T-1);
if ~exist(output_dir,'dir'), mkdir(output_dir); end
% Source terminal policies are calibrated pre-shock SS cpol/c_fine.
% Compare them BEFORE analyzing backwards recursion to distinguish boundary
% mismatches from EGM interpolation errors.
save_csv(output_dir,'fixed_src_terminal_coarse',s.cpol);
save_csv(output_dir,'fixed_src_terminal_fine',s.c_fine);
save_csv(output_dir,'fixed_src_rate',s.r_t_fixed(:));
save_csv(output_dir,'fixed_src_wage',s.w_nr_guess_fixed(:));
save_csv(output_dir,'fixed_src_profits',s.Profits_fixed(:));
save_csv(output_dir,'fixed_src_q',s.q_agg_fixed(:));
save_csv(output_dir,'fixed_src_initial_distribution',s.In_dis_fixed);
% Original asset residual is:
% ((A_household(t) - B(t))/q(t) - K(t+1)).
Knext=s.k_l_t_fixed(1:T-1);
assets=s.new_shock_b_agg(1:T-1) + ...
       s.q_agg_fixed(1:T-1).*(Knext+s.resid(1:T-1));
save_csv(output_dir,'fixed_src_household_assets',assets(:));
for t=[1 9 T-1]
    save_csv(output_dir,sprintf('fixed_src_consumption_t%d',t), ...
             s.c_tran_fine_fixed(:,:,t));
    save_csv(output_dir,sprintf('fixed_src_savings_t%d',t), ...
             s.dec_tran_fine_fixed(:,:,t));
end
fprintf('ORIGINAL FIXED AUDIT EXPORTED: saved rates, wages, profits, q, assets, distribution, date 1/9 fine policies\n');
fprintf('Original fixed source max raw 599 residual %.17g; no GE root-finding performed.\n', ...
        max(abs(s.resid(:))));
end

function save_csv(out,name,data)
dlmwrite(fullfile(out,[name '.csv']),data,'precision','%.17g');
end
