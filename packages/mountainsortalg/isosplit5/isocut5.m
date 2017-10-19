function [dip_score,cutpoint,info]=isocut5(samples,sample_weights,opts)
if nargin<1, test_isocut5; return; end;
if nargin<2, sample_weights=[]; end;
if (nargin<3), opts=struct; end;

if (~isfield(opts,'already_sorted')) opts.already_sorted=0; end;
if (~isfield(opts,'return_info')) opts.return_info=[]; end;
if (~isfield(opts,'num_bins_factor')) opts.num_bins_factor=1; end;

info=struct;

[~,N]=size(samples);
if (N==0)
    error('Error in isocut5: N is zero.');
end;
num_bins_factor=opts.num_bins_factor;
num_bins=ceil(sqrt(N/2)*num_bins_factor);

if (length(sample_weights)==0), sample_weights=ones(1,N); end;

if (opts.already_sorted)
    X=samples;
else
    [X,sort_inds]=sort(samples);
    sample_weights=sample_weights(sort_inds);
end;

while 1
    num_bins_1=ceil(num_bins/2);
    num_bins_2=num_bins-num_bins_1;
    intervals=[1:num_bins_1,num_bins_2:-1:1];
    alpha=(N-1)/sum(intervals);
    intervals=intervals*alpha;
    inds=floor([1,1+cumsum(intervals)]);
    N_sub=length(inds);
    if (min(intervals)>=1)
        break;
    else
        num_bins=num_bins-1;
    end;
end;

cumsum_sample_weights=cumsum(sample_weights);

X_sub=X(inds);
spacings=X_sub(2:end)-X_sub(1:end-1);
%multiplicities=inds(2:end)-inds(1:end-1);
multiplicities=cumsum_sample_weights(inds(2:end))-cumsum_sample_weights(inds(1:end-1));
densities=multiplicities./spacings;

densities_unimodal_fit=jisotonic5(densities,'updown',multiplicities);

[~,peak_density_ind]=max(densities_unimodal_fit); peak_density_ind=peak_density_ind(1);
[ks_left,ks_left_index]=compute_ks5(multiplicities(1:peak_density_ind),densities_unimodal_fit(1:peak_density_ind).*spacings(1:peak_density_ind));
[ks_right,ks_right_index]=compute_ks5(multiplicities(end:-1:peak_density_ind),densities_unimodal_fit(end:-1:peak_density_ind).*spacings(end:-1:peak_density_ind));
ks_right_index=length(spacings)-ks_right_index+1;

if (ks_left>ks_right)
    critical_range=1:ks_left_index;
    dip_score=ks_left;
else
    critical_range=ks_right_index:length(spacings);
    dip_score=ks_right;
end;

%ks=compute_ks4(multiplicities,densities_unimodal_fit./densities.*multiplicities);
%dip_score=ks;

densities_resid=densities-densities_unimodal_fit;
densities_resid_fit=jisotonic5(densities_resid(critical_range),'downup',spacings(critical_range));
[~,cutpoint_ind]=min(densities_resid_fit); cutpoint_ind=cutpoint_ind(1);
cutpoint_ind=critical_range(1)+cutpoint_ind-1;
cutpoint=(X_sub(cutpoint_ind)+X_sub(cutpoint_ind+1))/2;

if (opts.return_info)
    info.spacings=spacings;
    info.lefts=X_sub(1:end-1);
    info.rights=X_sub(2:end);
    info.centers=(info.lefts+info.rights)/2;
    info.densities=densities;
    info.densities_unimodal=densities_unimodal_fit;
    info.critical_range=critical_range;
    info.densities_bimodal=densities_resid_fit+densities_unimodal_fit(critical_range);
    info.plot_xx=zeros(1,(N_sub-1)*2);
    info.plot_xx(1:2:end)=info.lefts;
    info.plot_xx(2:2:end)=info.rights;
    info.plot_densities=zeros(1,(N_sub-1)*2);
    info.plot_densities(1:2:end)=info.densities;
    info.plot_densities(2:2:end)=info.densities;
    info.plot_densities_unimodal=zeros(1,(N_sub-1)*2);
    info.plot_densities_unimodal(1:2:end)=info.densities_unimodal;
    info.plot_densities_unimodal(2:2:end)=info.densities_unimodal;
    %info.plot_densities_bimodal=zeros(1,(N_sub-1)*2);
    %info.plot_densities_bimodal(1:2:end)=info.densities_bimodal;
    %info.plot_densities_bimodal(2:2:end)=info.densities_bimodal;
    [counts,info.hist_bins]=hist_with_weights(samples,sample_weights,num_bins*3);
    info.hist_counts=counts;
    info.hist_densities=counts/(info.hist_bins(2)-info.hist_bins(1));
else
    info=struct;
end;

function [ks]=compute_ks4(counts1,counts2)
S1=cumsum(counts1)/sum(counts1);
S2=cumsum(counts2)/sum(counts2);
ks=max(abs(S1-S2));
ks=ks*sqrt((sum(counts1)+sum(counts2))/2);

function [best_ks,best_len]=compute_ks5(counts1,counts2)
len=length(counts1);
best_ks=-inf;
while (len>=4)||(len==length(counts1))
    ks=compute_ks4(counts1(1:len),counts2(1:len));
    if (ks>best_ks)
        best_ks=ks;
        best_len=len;
    end;
    len=floor(len/2);
end;

function [counts,bins]=hist_with_weights(X,weights,num_bins)
bin_width=(max(X(:))-min(X(:)))/num_bins;
bin_ints=round(X/bin_width);
i1=min(bin_ints);
i2=max(bin_ints);
bins=(i1:i2)*bin_width;
ii=bin_ints-i1+1;
counts=accumarray(ii',weights',[length(bins),1]);

function test_isocut5
close all;

num_trials=100;
cutpoints=[];
dip_scores=[];
run_times=[];

compare_with_isocut5_mex=0;

cutpoints_mex=[];
dip_scores_mex=[];

for trial=1:num_trials;

N0=1e2;

X=[randn(1,10*N0),randn(1,3*N0)+5];
sample_weights=[1*ones(1,10*N0),1*ones(1,3*N0)];

opts.return_info=(trial==1);

tA=tic;
[dip_score,cutpoint,info]=isocut5(X,sample_weights,opts);
if compare_with_isocut5_mex
[dip_score2,cutpoint2]=isocut5_mex(X);
cutpoints_mex(end+1)=cutpoint2;
dip_scores_mex(end+1)=dip_score2;
end;
run_time=toc(tA);

cutpoints(end+1)=cutpoint;
dip_scores(end+1)=dip_score;
run_times(end+1)=run_time;

if (trial==1)
    X_first=X;
end;

if (trial==1)    
    figure; hold on;
    bar(info.hist_bins,info.hist_densities,'FaceColor',[0.8,0.8,0.8],'EdgeColor',[0.8,0.8,0.8]);
    plot(info.plot_xx,info.plot_densities,'k','LineWidth',2);
    plot(info.plot_xx,info.plot_densities_unimodal,'g','LineWidth',2);
    %plot(info.plot_xx,info.plot_densities_bimodal,'b','LineWidth',2);
    vline0(cutpoint);
    title(sprintf('dip score = %g, N=%g',dip_score,length(X)));
    drawnow;
end;

end;

figure;
subplot(1,3,1);
hist(dip_scores,num_trials);
xlim([min(0,min(dip_scores)),max(dip_scores)]);
title(sprintf('dip scores - avg = %g',mean(dip_scores)));
subplot(1,3,2);
hist(cutpoints,num_trials);
xlim([min(X_first),max(X_first)]);
title(sprintf('cutpoints - avg = %g',mean(cutpoints)));
subplot(1,3,3);
hist(run_times,num_trials);
title(sprintf('run times - avg = %g seconds',round(mean(run_times)*100000)/100000));
set(gcf,'position',[50,50,1000,500]);

if compare_with_isocut5_mex
figure;
plot(cutpoints,cutpoints_mex,'b.');
title('cutpoint comparison');
figure;
plot(dip_scores,dip_scores_mex,'b.');
title('dip scores comparison');
end;
