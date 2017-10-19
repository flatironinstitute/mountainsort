function [labels,info]=isosplit5(X,opts)
% isosplit5 - perform clustering using isotonic regression (jfm, may 2015 - dec 2016)
%
% labels = isosplit5(X,opts) 
%   X is M x N, M=#dimensions, N=#samples
%   labels: 1xN vector of labels from 1..K,  K = # clusters
%
%   opts.isocut_threshold -- threshold for determining isocut tests
%   opts.min_cluster_size -- minimum cluster size
%   opts.K_init -- number of clusters in initial parcelation
%   opts.refine_clusters -- whether to recursively apply isosplit to refine
%        clusters
%   opts.max_iterations -- maximum number of iterations before stopping
%   opts.verbose
%   opts.verbose_pause_duration
%   opts.whiten_cluster_pairs -- whether to whiten at each comparison
%   opts.initial_labels -- optional -- if provided will skip the parcelation step
%
% Magland 5/19/2015, updated dec 2016

%% Default parameters and self test
if nargin<1, test_isosplit5; return; end;
if nargin<2, opts=struct; end;
if ~isfield(opts,'isocut_threshold') opts.isocut_threshold=1; end;
if ~isfield(opts,'min_cluster_size') opts.min_cluster_size=10; end;
if ~isfield(opts,'K_init') opts.K_init=200; end;
if ~isfield(opts,'refine_clusters'), opts.refine_clusters=true; end;
if ~isfield(opts,'max_iterations_per_pass'), opts.max_iterations_per_pass=500; end;
if ~isfield(opts,'verbose') opts.verbose=0; end;
if ~isfield(opts,'verbose_isocut') opts.verbose_isocut=0; end;
if ~isfield(opts,'verbose_pause_duration') opts.verbose_pause_duration=0.5; end;
if ~isfield(opts,'whiten_cluster_pairs') opts.whiten_cluster_pairs=1; end;
if ~isfield(opts,'initial_labels') opts.initial_labels=[]; end;
if ~isfield(opts,'prevent_merge') opts.prevent_merge=false; end;
if ~isfield(opts,'one_comparison_at_a_time') opts.one_comparison_at_a_time=false; end; % for verbose=1
if ~isfield(opts,'return_iterations') opts.return_iterations=0; end;

%% Initialize the timers for diagnostic
timers.get_pairs_to_compare=0;
timers.compare_pairs=0;
timers.compute_centers=0;

[M,N]=size(X);
info=struct;
if (opts.return_iterations)
    info.iterations={};
end;

%% Compute the initial clusters
if (isempty(opts.initial_labels))
    target_parcel_size=opts.min_cluster_size;
    target_num_parcels=opts.K_init;
    % !! important not to do a final reassign because then the shapes will not
    % be conducive to isosplit iterations -- hexagons are not good for isosplit!
    data.labels=parcelate2(X,target_parcel_size,target_num_parcels,struct('final_reassign',0));
    Kmax=max(data.labels);
else
    data.labels=opts.initial_labels;
    Kmax=max(data.labels);
    inds0=find(data.labels==0);
    if (~isempty(inds0))
        inds1=find(data.labels>0);
        if (~isempty(inds1))
            idx1=knnsearch(X(:,inds1)',X(:,inds0)');
            data.labels(inds0)=data.labels(inds1(idx1));
        else
            data.labels(inds0)=1;
        end;
    end;
end;


%debug
%labels=data.labels;
%return;

%% Compute the cluster centers
ttt=tic;
data.centers=compute_centers(X,data.labels);
timers.compute_centers=timers.compute_centers+toc(ttt);

if (opts.verbose)
    fA=figure;
    info.frames={};
end;

%% Repeat while something has been merged in the pass
final_pass=false; % plus we do one final pass at the end
data.comparisons_made=zeros(Kmax,Kmax); % Keep a matrix of comparisons that have been made in this pass
while 1 % Passes
    something_merged=false; % Keep track of whether something has merged in this pass. If not, do a final pass.
    
    clusters_changed_vec_in_pass=zeros(1,Kmax); % Keep track of the clusters that have changed in this pass so that we can update the comparisons_made matrix at the end
    iteration_number=0;
    while 1 % Iterations
        iteration_number=iteration_number+1;
        if (iteration_number>opts.max_iterations_per_pass)
            error('max iterations per pass exceeded');
            break;
        end;
        % The active labels are those that are still being used
        active_labels_vec=zeros(1,Kmax);
        active_labels_vec(data.labels)=1;
        active_labels=find(active_labels_vec);
        active_centers=data.centers(:,active_labels);

        % Find the pairs to compare on this iteration
        % These will be closest pairs of active clusters that have not yet
        % been compared in this pass
        ttt=tic;
        [inds1,inds2]=get_pairs_to_compare(active_centers,data.comparisons_made(active_labels,active_labels),opts);
        timers.get_pairs_to_compare=timers.get_pairs_to_compare+toc(ttt);

        % If we didn't find any, break from this iteration
        if (length(inds1)==0)
            % Nothing else to compare.
            break;
        end;
        
        % Show the clusters if we are in verbose mode
        if (opts.verbose)
%             labels_map=zeros(1,Kmax);
%             active_labels_vec=zeros(1,Kmax);
%             active_labels_vec(data.labels)=1;
%             active_labels=find(active_labels_vec);
%             for ii=1:length(active_labels)
%                 labels_map(active_labels(ii))=ii;
%             end;
%             labels_mapped=labels_map(data.labels);
            ooo.draw_axes=false;
            ooo.draw_legend=false;
            if (length(unique(data.labels))<=10)
                ooo.draw_legend=true;
            end;
            figure(fA); ms_view_clusters_0(X(:,:),data.labels,ooo);
            %title(sprintf('iteration %d',iteration_number));
            set(gca,'xtick',[]); set(gca,'ytick',[]);
            pause(opts.verbose_pause_duration);
            info.frames{end+1}=getframe;
        end;
        
        old_labels=data.labels; % So we can determine the number of label changes for diagnostics

        % Actually compare the pairs -- in principle this operation could be parallelized
        ttt=tic;
        [data.labels,clusters_changed,compare_pairs_info]=compare_pairs(X,data.labels,active_labels(inds1),active_labels(inds2),opts);
        clusters_changed_vec_in_pass(clusters_changed)=1;
        timers.compare_pairs=timers.compare_pairs+toc(ttt);
        
        if (opts.return_iterations)
            AA.pairs_to_compare_1=active_labels(inds1);
            AA.pairs_to_compare_2=active_labels(inds2);
            AA.old_labels=old_labels;
            AA.new_labels=data.labels;
            AA.compare_pairs_info=compare_pairs_info;
            info.iterations{end+1}=AA;
        end;

        % Update which comparisons have been made
        for j=1:length(inds1)
            data.comparisons_made(active_labels(inds1(j)),active_labels(inds2(j)))=1;
            data.comparisons_made(active_labels(inds2(j)),active_labels(inds1(j)))=1;
        end;

        % Recompute the centers -- note: maybe this should only apply to those that changed?
        ttt=tic;
        data.centers=compute_centers(X,data.labels);
        timers.compute_centers=timers.compute_centers+toc(ttt);

        % For diagnostics, cound the number of changes
        total_num_label_changes=length(find(data.labels~=old_labels));
    
        % Determine whether something has merged
        new_active_labels_vec=zeros(1,N);
        new_active_labels_vec(data.labels)=1;
        new_active_labels=find(new_active_labels_vec);
        if (length(new_active_labels)<length(active_labels))
            something_merged=1;
        end;
        
        %break;
    end;
    % zero out the comparisons made matrix only for those that have changed
    clusters_changed=find(clusters_changed_vec_in_pass);
    for j=1:length(clusters_changed)
        data.comparisons_made(clusters_changed(j),:)=0;
        data.comparisons_made(:,clusters_changed(j))=0;
    end;
    
    if (something_merged) final_pass=false; end;
    if (final_pass) break; end; % This was the final pass and nothing has merged
    if (~something_merged) final_pass=true; end; % If we are done, do one last pass for final redistributes
    
end;

% This is the result
labels=data.labels;

% But we should remap the labels to occupy the first natural numbers
labels_map=zeros(1,Kmax);
active_labels_vec=zeros(1,Kmax);
active_labels_vec(labels)=1;
active_labels=find(active_labels_vec);
for ii=1:length(active_labels)
    labels_map(active_labels(ii))=ii;
end;
labels=labels_map(labels);

% Return the timers in the info
info.timers=timers;

% If the user wants to refine the clusters, then we repeat isosplit on each
% of the new clusters, recursively. Unless we only found one cluster.
if ((opts.refine_clusters)&&(max(labels)>1))
    opts2=opts;
    opts2.initial_labels=[];
    opts2.refine_clusters=true; % Maybe we should provide an option on whether to do recursive refinement
    K=max(labels);
    labels_split=zeros(1,N);
    for k=1:K
        inds_k=find(labels==k);
        X_k=X(:,inds_k);
        labels_k=isosplit5(X_k,opts2);
        labels_split(inds_k)=max(labels_split)+labels_k;
    end;
    labels=labels_split;
end;

function centers=compute_centers(X,labels)
[M,N]=size(X);
centers=zeros(M,N);
counts=accumarray(labels',1,[N,1])';
for m=1:M
    centers(m,:)=accumarray(labels',X(m,:)',[N,1])';
end;
centers(:,find(counts))=centers(:,find(counts))./repmat(counts(find(counts)),M,1);

function [new_labels,clusters_changed,info]=compare_pairs(X,labels,k1s,k2s,opts)
info=struct;
if (opts.return_iterations)
    info.merge_tests={};
end;
clusters_changed_vec=zeros(1,max(labels));
new_labels=labels;
for i1=1:length(k1s)
    k1=k1s(i1);
    k2=k2s(i1);
    inds1=find(labels==k1);
    inds2=find(labels==k2);
    if ((length(inds1)>0)&&(length(inds2)>0))
        if ((length(inds1)<opts.min_cluster_size)||(length(inds2)<opts.min_cluster_size))
            do_merge=1;
        else
            inds12=cat(2,inds1,inds2);
            L12_old=cat(2,ones(1,length(inds1)),2*ones(1,length(inds2)));
            [do_merge,L12,proj,cutpoint]=merge_test(X(:,inds1),X(:,inds2),opts);
            if (opts.return_iterations)
                MM.do_merge=do_merge;
                MM.L12=L12;
                MM.proj=proj;
                MM.cutpoint=cutpoint;
                MM.k1=k1;
                MM.k2=k2;
                info.merge_tests{end+1}=MM;
            end;
            
            if (opts.verbose_isocut)
                title(sprintf('Compare %d %d',k1,k2));
            end;
        end;
        if (do_merge)
            if (~opts.prevent_merge)
                new_labels(find(new_labels==k2))=k1;
                clusters_changed_vec(k1)=1;
                clusters_changed_vec(k2)=1;
            end;
        else
            %redistribute
            new_labels(inds12(find(L12==1)))=k1;
            new_labels(inds12(find(L12==2)))=k2;
            if (length(find(L12~=L12_old))>0)
                clusters_changed_vec(k1)=1;
                clusters_changed_vec(k2)=1;
            end;
        end;
    end;
end;
clusters_changed=find(clusters_changed_vec);

function [X1b,X2b,V]=whiten_two_clusters_b(X1,X2)
M=size(X1,1);
N1=size(X1,2);
N2=size(X2,2);

centroid1=mean(X1,2);
centroid2=mean(X2,2);

X1_centered=X1-repmat(centroid1,1,N1);
X2_centered=X2-repmat(centroid2,1,N2);

C1=(X1_centered*X1_centered')/N1;
C2=(X2_centered*X2_centered')/N2;
avg_cov=(C1+C2)/2;

X1b=X1;
X2b=X2;
V=centroid2-centroid1;
if (abs(det(avg_cov))>1e-6)
    inv_avg_cov=inv(avg_cov);
    V=inv_avg_cov*V;
end;
V=V/sqrt(V'*V);

function [X1b,X2b,V]=whiten_two_clusters(X1,X2)
M=size(X1,1);
N1=size(X1,2);
N2=size(X2,2);

% Important to subtract the two centroids before whitening!
centroid1=mean(X1,2);
centroid2=mean(X2,2);
Y1=X1-repmat(centroid1,1,N1);
Y2=X2-repmat(centroid2,1,N2);

% Combine the data
Y=cat(2,Y1,Y2);
N=N1+N2;

% Obtain the whitening matrix using svd
if (N>=M)
    [U,D,V] = svd(Y,'econ');
    D(D~=0)=1./D(D~=0);
    % Amd apply it to the original (non-mean subtracted) data
    X1b=sqrt(N-1)*U*D(1:M,1:M)*(U'*X1);
    X2b=sqrt(N-1)*U*D(1:M,1:M)*(U'*X2);
else
    %too few points to whiten
    X1b=X1;
    X2b=X2;
end;

% The best direction is now the one connecting the centroids.
centroid1b=mean(X1b,2);
centroid2b=mean(X2b,2);
V=centroid2b-centroid1b;

function [do_merge,new_labels,projection12,cutpoint]=merge_test(X1_in,X2_in,opts)
if opts.whiten_cluster_pairs
    [X1,X2,V]=whiten_two_clusters_b(X1_in,X2_in);
else
    X1=X1_in;
    X2=X2_in;
    centroid1=mean(X1,2);
    centroid2=mean(X2,2);
    V=centroid2-centroid1;
    V=V/sqrt(V'*V);
end;

[~,N1]=size(X1); [~,N2]=size(X2);
if ((N1==0)||(N2==0))
    error('Error in merge test: N1 or N2 is zero');
end;

projection1=V'*X1;
projection2=V'*X2;
projection12=cat(2,projection1,projection2);
[dipscore,cutpoint]=isocut5(projection12,ones(size(projection12)));
do_merge=(dipscore<opts.isocut_threshold);
if (opts.verbose_isocut)
    figure; view_isocut_hist(projection12,do_merge,cutpoint);
end;
%cutpoint=isocut(projection12,opts.isocut_threshold);
%ret=(cutpoint~=0);
new_labels=ones(1,N1+N2);
new_labels(find(projection12>=cutpoint))=2;

function dists=make_dists_matrix(centers)
[M,N]=size(centers);
dists=zeros(N,N);
[aa,bb]=ndgrid(1:N,1:N);
for m=1:M
    dists=dists+reshape((centers(m,aa(:))-centers(m,bb(:))).^2,N,N);
end;
dists=sqrt(dists);

function [inds1,inds2]=get_pairs_to_compare(centers,comparisons_made,opts)
[M,N]=size(centers);
inds1=[];
inds2=[];
pair_dists=[];
dists=make_dists_matrix(centers);
dists(find(comparisons_made(:)))=inf;
for j=1:N
    dists(j,j)=inf;
end;
% important to only take the mutal closest pairs -- unlike how we originally did it
%something_changed=1;
%while (something_changed)
    %something_changed=0;
    [~,best_inds]=min(dists,[],1);
    for j=1:N
        if (best_inds(j)>j)
            if (best_inds(best_inds(j))==j) % mutual
                if (dists(j,best_inds(j))<inf)
                    inds1(end+1)=j;
                    inds2(end+1)=best_inds(j);
                    pair_dists(end+1)=dists(j,best_inds(j));
                    dists(j,:)=inf;
                    dists(:,j)=inf;
                    dists(best_inds(j),:)=inf;
                    dists(:,best_inds(j))=inf;
                    %something_changed=1;
                end;
            end;
        end;        
    end;
%end;

if (opts.one_comparison_at_a_time)
    if (length(inds1)>0)
        [~,ii]=min(pair_dists); ii=ii(1);
        inds1=inds1(ii);
        inds2=inds2(ii);
    end;
end;

function view_isocut_hist(X,do_merge,cutpoint)
if (do_merge)
    inds1=1:length(X);
    inds2=[];
else
    inds1=find(X<cutpoint);
    inds2=find(X>=cutpoint);
end;
[counts,bins]=hist(X,ceil(length(X)/10));
counts1=hist(X(inds1),bins);
counts2=hist(X(inds2),bins);
bar(bins,counts,1,'FaceColor',[0.8,0.8,0.8],'EdgeColor',[0.8,0.8,0.8]);
hold on;
bar(bins,counts1,1,'FaceColor',[0.3,0.3,0.3],'EdgeColor',[0.3,0.3,0.3]);
bar(bins,counts2,1,'FaceColor',[0.6,0.6,0.6],'EdgeColor',[0.6,0.6,0.6]);
hold off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,labels]=create_multimodal_nd(A)
M=length(A{1}.center);
X=zeros(M,0);
labels=zeros(1,0);
for j=1:length(A)
    A0=A{j};
    tmp=randn(M,A0.N);
    %b=[A0(4),A0(6);-A0(6),A0(5)];
    tmp=A0.cov*tmp;
    for m=1:M
        tmp(m,:)=A0.center(m)+tmp(m,:);
    end;
    X=cat(2,X,tmp);
    labels=cat(2,labels,j*ones(1,A0.N));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,labels]=generate_dataset

N0=1e4;
num_noise_dims=0;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[6,0]; A2.cov=[1,0;0,2.5];
%A3.N=N0/2; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
%A4.N=ceil(N0/16); A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
A3.N=1000; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
A4.N=500; A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];

A5.N=N0/2; A5.center=[11,7]; A5.cov=[1,0;0,2.5];

AA={A1,A2,A3,A4,A5};
for j=1:length(AA)
    M=length(AA{j}.center);
    center2=rand(1,M+num_noise_dims)*0;
    center2(1:M)=AA{j}.center;
    AA{j}.center=center2;
    cov2=eye(M+num_noise_dims);
    cov2(1:M,1:M)=AA{j}.cov;
    AA{j}.cov=cov2;
end;
[X,labels]=create_multimodal_nd(AA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_isosplit5
rng(104);
close all;

[X,true_labels]=generate_dataset;
figure; ms_view_clusters_0(X(1:2,:),true_labels);
title('Truth');

ttt=tic;
[labels2,info]=isosplit5(X,struct('verbose',1,'refine_clusters',0,'whiten_cluster_pairs',1));
fprintf('Time for isosplit5: %g\n',toc(ttt));
figure; ms_view_clusters_0(X(1:2,:),labels2);
title('isosplit5');
disp(info.timers);

% try
%     test1=isosplit5_mex(randn(2,100));
% catch
%     compile_mex_isosplit5;
%     compile_mex_isocut5;
% end

ttt=tic;
labels_mex=isosplit5_mex(X);
fprintf('Time for isosplit5_mex: %g\n',toc(ttt));
figure; ms_view_clusters_0(X(1:2,:),labels_mex);
title('isosplit5 mex');
