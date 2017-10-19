function [labels,parcels]=parcelate2(X,target_parcel_size,target_num_parcels,opts)
if nargin<1, test_parcelate2; return; end;
if nargin<4, opts=struct; end;
% final_reassign requires knnsearch - which requires a toolbox
if ~isfield(opts,'final_reassign'), opts.final_reassign=0; end;

[M,N]=size(X);
labels=zeros(1,N);

parcels={};

% Begin with a single parcel with all the points
P.indices=1:N;
P.centroid=mean(X(:,P.indices),2);
P.radius=compute_max_distance(P.centroid,X(:,P.indices));
labels(P.indices)=1;
parcels{end+1}=P;

%fprintf('Radius of initial cluster is %g\n',P.radius);

split_factor=3; % split factor around 2.71 is in a sense ideal

% continue dividing while we have fewer than the target number of parcels
while (length(parcels)<target_num_parcels)
    parcel_sizes=get_parcel_sizes(parcels);
    parcel_radii=get_parcel_radii(parcels);
    if (length(find((parcel_radii>0)&(parcel_sizes>target_parcel_size)))==0)
        % nothing else will ever be split
        break;
    end;
    % the target radius is a bit smaller than the largest parcel radius
    target_radius=max(parcel_radii(find(parcel_sizes>target_parcel_size)))*0.95;
    
    % step through each parcel
    something_changed=0;
    p_index=1;
    while (p_index<=length(parcels))
        inds=parcels{p_index}.indices;
        rad=parcels{p_index}.radius;
        sz=length(inds);
        if (sz>target_parcel_size)&&(rad>=target_radius)
            % the parcel has more than the tnarget number of points and is larger than the target radius
            %iii=randsample(sz,split_factor);
            iii=1:split_factor;
            pts=X(:,inds(iii)); % take the first few points as the parcel representatives (deterministic)
            dm=distance_matrix(pts,X(:,inds));
            [~,assignments]=min(dm,[],1); % make assignments based on distances to these points
            if ((length(find(assignments==1))>0)&&(length(find(assignments~=1))>0)) % make sure something but not everything got assigned to 1                
                something_changed=1;
                parcels{p_index}.indices=inds(find(assignments==1));
                parcels{p_index}.centroid=mean(X(:,parcels{p_index}.indices),2);
                parcels{p_index}.radius=compute_max_distance(parcels{p_index}.centroid,X(:,parcels{p_index}.indices));
                labels(parcels{p_index}.indices)=p_index;
                for jj=2:split_factor
                    PP.indices=inds(find(assignments==jj));
                    PP.centroid=mean(X(:,PP.indices),2);
                    PP.radius=compute_max_distance(PP.centroid,X(:,PP.indices));
                    parcels{end+1}=PP;
                    labels(PP.indices)=length(parcels);
                end;
            else
                warning('There was an issue splitting a parcel. This could result from the same point being included multiple times in the data.');
            end;
        else
            p_index=p_index+1;
        end;
    end;
    if (~something_changed)
        break;
    end;
end;

if (opts.final_reassign)
    centroids=get_parcel_centroids(parcels);
    labels=knnsearch(centroids',X','K',1)';
end;
%dm=distance_matrix(centroids,X);
%[~,labels]=min(dm,[],1);

function ret=get_parcel_sizes(parcels)
ret=zeros(1,length(parcels));
for j=1:length(parcels)
    ret(j)=length(parcels{j}.indices);
end;

function ret=get_parcel_radii(parcels)
ret=zeros(1,length(parcels));
for j=1:length(parcels)
    ret(j)=parcels{j}.radius;
end;

function ret=get_parcel_centroids(parcels)
ret=zeros(size(parcels{1}.centroid,1),length(parcels));
for j=1:length(parcels)
    ret(:,j)=parcels{j}.centroid;
end;

function ret=compute_max_distance(pt,X)
dm=distance_matrix(pt,X);
ret=max(dm);
if (length(ret)==0)
    ret=0;
end;

function ret=distance_matrix(X1,X2)
[M,N1]=size(X1);
[M,N2]=size(X2);
ret=zeros(N1,N2);
[aa,bb]=ndgrid(1:N1,1:N2);
for m=1:M
    ret=ret+reshape((X1(m,aa(:))-X2(m,bb(:))).^2,N1,N2);
end;
ret=sqrt(ret);

function test_parcelate2

target_parcel_size=100;
target_num_parcels=500;
M=2;
N=1e6;
X=randn(M,N);
tic;
labels=parcelate2(X,target_parcel_size,target_num_parcels);
toc
figure; ms_view_clusters(X,labels);

