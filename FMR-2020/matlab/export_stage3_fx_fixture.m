function export_stage3_fx_fixture(original_root,input_dir,output_dir)
% Original source portfolio revaluation and basefun distribution splitting,
% using exactly the same pre-shock household distribution as Julia.
if nargin~=3, error('Provide original_root, input_dir, output_dir'); end
addpath(fullfile(original_root,'Exogenous_files'));
if ~exist('basefun','file'), error('Original basefun.m not found'); end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
grid=dlmread(fullfile(input_dir,'fine_grid.csv'),',');
G0=dlmread(fullfile(input_dir,'mu_pre.csv'),',');
meta=dlmread(fullfile(input_dir,'metadata.csv'),',');
n=length(grid);
assert(size(G0,1)==n && size(G0,2)==7);
assert(size(meta,1)==2 && size(meta,2)==11);
roots=zeros(2,1);
for case_id=1:2
    dat=meta(case_id,:);
    K=dat(1); qprev=dat(2); Bprev=dat(3);
    xprev=dat(4); hat=dat(5);
    C=dat(6); E=dat(7); D=dat(8); theta=dat(9);
    author_eq=@(x) C+E/x-D/x^theta;
    [result,err,exitflag]=fsolve(author_eq,0.8,optimset('Display','off'));
    if ~(exitflag>0 && result>0 && abs(err)<1e-8)
        error('Author original scalar trade fsolve did not converge');
    end
    roots(case_id)=result;
    netrich=sum(sum(G0(grid>0,:),2).*grid(grid>0));
    leverage=(1-hat)*K*qprev/netrich;
    fixed_capital=hat*K;
    for branch=1:2
        xnew=dat(9+branch);
        Bgrid=-(max(leverage*grid+fixed_capital,0)-grid);
        Kgrid=max(leverage*grid+fixed_capital,0);
        newassets=Kgrid+Bgrid/(xnew/xprev);
        G=zeros(size(G0));
        for ai=1:n
            [vals,inds]=basefun(grid,n,newassets(ai));
            for j=1:size(G0,2)
                G(inds(1),j)=G(inds(1),j)+G0(ai,j)*vals(1);
                G(inds(2),j)=G(inds(2),j)+G0(ai,j)*vals(2);
            end
        end
        if case_id==1, label='unleveraged'; else, label='leveraged'; end
        if branch==1, rootlabel='low'; else, rootlabel='high'; end
        dlmwrite(fullfile(output_dir,[label '_' rootlabel '_distribution.csv']), ...
                 G,'precision','%.17g');
        fprintf('SOURCE FX leverage=%g root=%s mass=%g wealth=%g\n',...
                hat,rootlabel,sum(G(:)),sum(G'*grid));
    end
end
dlmwrite(fullfile(output_dir,'original_fsolve_root.csv'),roots,...
         'precision','%.17g');
fprintf('ORIGINAL FX SOURCE-BASED INTERPOLATION COMPLETE\n');
end
