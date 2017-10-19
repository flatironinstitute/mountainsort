function [B,MSEs]=jisotonic5(A,direction,weights)
% jisotonic5 - isotonic regression (jfm, may 2105)
%
% [B,MSEs] = jisotonic5(A,direction,weights)
%   A is the input a vector
%   direction = 'increading', 'decreasing', 'updown', or 'downup'
%   weights is the optional input weight vector (same size as A)
%   B is the output vector (same size as A)
%   MSEs is used internally for 'updown' and 'downup' directions
%
% Magland 5/19/2015

if (nargin<1)
	jisotonic5_test; %run the test code
	return;
end;
if (nargin<2)
	direction='increasing';
end;
if (~isrow(A)) A=A'; end;
if (nargin<3)
	weights=ones(size(A));
end;

try
    [test1,test2]=jisotonic5_mex(1,1);
catch
    compile_mex_jisotonic5;
end;

if (strcmp(direction,'decreasing'))
	%[B,MSEs]=jisotonic5(-A,'increasing',weights); B=-B;
    [B,MSEs]=jisotonic5_mex(-A,weights); B=-B;
	return;
elseif (strcmp(direction,'updown'))
	%[B1,MSE1]=jisotonic5(A,'increasing',weights);
    [B1,MSE1]=jisotonic5_mex(A,weights);
	%[B2,MSE2]=jisotonic5(A(end:-1:1),'increasing',weights(end:-1:1));
    [B2,MSE2]=jisotonic5_mex(A(end:-1:1),weights(end:-1:1));
	B2=B2(end:-1:1);
	MSE2=MSE2(end:-1:1);
	MSE0=MSE1+MSE2;
	
	[~,best_ind]=min(MSE0);
	%C1=jisotonic5(A(1:best_ind),'increasing',weights(1:best_ind));
    [C1,dum]=jisotonic5_mex(A(1:best_ind),weights(1:best_ind));
	%C2=jisotonic5(A(best_ind:end),'decreasing',weights(best_ind:end));
    [C2,dum]=jisotonic5_mex(-A(best_ind:end),weights(best_ind:end)); C2=-C2;
	B=[C1(1:best_ind),C2(2:end)];
	if (isnan(B(1)))
		warning('jisotonic5: NaN');
	end;
	return;
elseif (strcmp(direction,'downup'))
	B=-jisotonic5(-A,'updown',weights);
	return;
else
	if (~strcmp(direction,'increasing'))
		error(['invalid direction in jisotonic5: ',direction]);
		return;
	end;
end;

try
[B,MSEs]=jisotonic5_mex(A,weights);
catch
	error('Unable to run mex file -- You must compile using: mex jisotonic5_mex.cpp');
	%mex(sprintf('%s/jisotonic5_mex.cpp',fileparts(mfilename('fullpath'))));
	%[B,MSEs]=jisotonic5_mex(A,weights);
end;

%assume increasing
% X=cell(1,0);
% N=length(A);
% 
% tmp.unweightedcount=1;
% tmp.count=weights(1);
% tmp.sum=A(1)*weights(1);
% tmp.sumsqr=A(1)^2*weights(1);
% lastind=1;
% X{lastind}=tmp;
% 
% MSEs=zeros(1,N);
% MSEs(1)=0;
% 
% for j=2:N
% 	tmp.unweightedcount=1;
% 	tmp.count=weights(j);
% 	tmp.sum=A(j)*weights(j);
% 	tmp.sumsqr=A(j)^2*weights(j);
% 	X{lastind+1}=tmp;
% 	MSEs(j)=MSEs(j-1);
% 	lastind=lastind+1;
% 	
% 	while true
% 		if (lastind<=1) break; end;
% 		tmp1=X{lastind-1};
% 		tmp2=X{lastind};
% 		prevMSE=tmp1.sumsqr-tmp1.sum^2/tmp1.count + tmp2.sumsqr-tmp2.sum^2/tmp2.count;
% 		if (tmp1.sum/tmp1.count<tmp2.sum/tmp2.count)
% 			break;
% 		else
% 			tmp.unweightedcount=tmp1.unweightedcount+tmp2.unweightedcount;
% 			tmp.count=tmp1.count+tmp2.count;
% 			tmp.sum=tmp1.sum+tmp2.sum;
% 			tmp.sumsqr=tmp1.sumsqr+tmp2.sumsqr;
% 			X{lastind-1}=tmp;
% 			lastind=lastind-1;
% 			newMSE=tmp.sumsqr-tmp.sum^2/tmp.count;
% 			MSEs(j)=MSEs(j)+newMSE-prevMSE;
% 		end;
% 	end;
% end;
% 
% B=zeros(1,N);
% ii=1;
% for k=1:lastind
% 	tmp0=X{k};
% 	for cc=1:tmp0.unweightedcount;
% 		B(ii+cc-1)=tmp0.sum/tmp0.count;
% 	end;
% 	ii=ii+tmp0.unweightedcount;
% end;

end

function jisotonic5_test

A=[1.1,3.2,5.3,2.4,4.5,6.6,3.7,5.8,7.9,10.0,9.2,11.2,7.3,3.4,5.5,2.6,3.7,1.8];
[B,BMSE]=jisotonic5(A);
[C,CMSE]=jisotonic5(A,'decreasing');
D=jisotonic5(A,'updown');
figure; plot(1:length(A),A,'b',1:length(B),B,'r',1:length(C),C,'g',1:length(D),D,'k');
figure; plot(1:length(BMSE),BMSE,'k');

end
