clear
clc

vars0=importdata('MF_psy_vars.mat');
NET=importdata('NET.mat');
load('SMvars.mat')
BMvars=importdata('SDR_NET1_PCnew.mat');

Nkeep=100;

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

% Remove bad variables
badvars=[];
for i=1:size(SMvars,2)
  Y=SMvars(:,i); grotKEEP=~isnan(Y);  
  %change the creterion when # of subjects changes
  if (sum(grotKEEP)>400) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.99)
      % the 3rd thing above is:  is the size of the largest equal-values-group too large?
    i=i; % do nothing
  else
    [i sum(grotKEEP)  max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))]
    badvars=[badvars i];
  end 
end

varskeep=setdiff([1:size(SMvars,2)], badvars);  
SMvars=SMvars(:,varskeep);

% Normalise and deconfound the original SM
varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grot=varsgrot; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
varsdCOV = (grot*grot') ./ (grotI*grotI');
%varsdCOV=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[uu2,dd]=eigs(varsdCOV,55);  % SVD (eigs actually)
SMpc=uu2*sqrt(dd)-conf*(pinv(conf)*(uu2*sqrt(dd)));    % deconfound again just to be safe

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu,ss,vv]=nets_svds(NETd,100); % SVD reduction
BMpc=uu*ss;

% BMvars=nets_demean(BMvars-conf*(pinv(conf)*BMvars));   

% %Apply PCA to further reduce the dimension of SDR BM
% [uu,ss,vv]=nets_svds(BMvars,Nkeep); % SVD reduction
% BMpc=uu*ss;
 
%%% CCA/
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(BMpc,SMpc);


%%% amount of variance between explained by canonical variates in original variables
for j=1:5
    TVE_SM(j,:)=corr(grotV(:,j),varsgrot,'rows','pairwise').^2;%varsgrot
    TVE_Net(j,:)=corr(grotU(:,j),NETd).^2; %NETd(:,1:size(NET,2))
end
AVE_SM=mean(TVE_SM,2,'omitnan');
AVE_Net=mean(TVE_Net,2,'omitnan');

for i=1:length(AVE_SM)
    cum(i)=sum(AVE_SM(1:i)');  
end


varsgrot(isnan(varsgrot))=0;
% Alternatively
b1=grotV(:,1:3)\varsgrot;
yhat=grotV(:,1:3)*b1;
Rsq = 1 - sum((varsgrot - yhat).^2)/sum((varsgrot).^2)

b2=grotU(:,1:3)\NETd;
yhat2=grotU(:,1:3)*b2;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)


