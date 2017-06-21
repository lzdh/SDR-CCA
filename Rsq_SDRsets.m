clear
clc

smpca100=importdata('pca100_smdata.mat');
smpca55=importdata('pca55_smdata.mat');
vars0=importdata('MF_psy_vars.mat');
vars=importdata('SDR_vars_PC.mat');
NET=importdata('SDR_NET1_PCnew.mat');
NET0=importdata('NET.mat');
load('SMvars.mat')
SM_header=importdata('SM_header.txt');

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3)]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

badvars=[];
for i=1:size(SMvars,2)
  Y=SMvars(:,i); grotKEEP=~isnan(Y);  
  %change the creterion when # of subjects changes
  if (sum(grotKEEP)>400) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95)
      % the 3rd thing above is:  is the size of the largest equal-values-group too large?
    i=i; % do nothing
  else
    [i sum(grotKEEP)  max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))]
    badvars=[badvars i];
  end 
end

varskeep=setdiff([1:size(SMvars,2)], badvars);  
SMvars=SMvars(:,varskeep);
SM_header=SM_header(varskeep');

varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end

vars=palm_inormal(vars);
for i=1:size(vars,2)
  grot=(isnan(vars(:,i))==0); grotconf=nets_demean(conf(grot,:)); vars(grot,i)=nets_normalise(vars(grot,i)-grotconf*(pinv(grotconf)*vars(grot,i)));
end

NET1=nets_demean(NET0);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
NET=nets_demean(NET-conf*(pinv(conf)*NET)); 

[uu1,ss1,vv1]=nets_svds(NETd,100); % SVD reduction
bmpca100=uu1*ss1;

[uu1,ss1,vv1]=nets_svds(NETd,1571); % SVD reduction
bmpca1571=uu1*ss1;


%Apply PCA to further reduce the dimension of SDR BM
[uu,ss,vv]=nets_svds(NET,100); % SVD reduction
BMpc=uu*ss;

vars(isnan(vars))=0;
varsgrot(isnan(varsgrot))=0;

b1=vars\varsgrot;
yhat=vars*b1;
Rsq = 1 - sum((varsgrot - yhat).^2)/sum((varsgrot).^2)

bsm_pca100=smpca100\varsgrot;
ypca100=smpca100*bsm_pca100;
Rsq_pca100 = 1 - sum((varsgrot - ypca100).^2)/sum((varsgrot).^2)

bsm_pca55=smpca55\varsgrot;
ypca55=smpca55*bsm_pca55;
Rsq_pca55 = 1 - sum((varsgrot - ypca55).^2)/sum((varsgrot).^2)


bbm=NET\NETd;
yhat2=NET*bbm;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)

bbm_sdrpca100=BMpc\NETd;
ybmsdrpca100=BMpc*bbm_sdrpca100;
Rsq2_sdrpca100 = 1 - sum((NETd - ybmsdrpca100).^2)/sum((NETd).^2)

bbm_pca100=bmpca100\NETd;
ybmpca100=bmpca100*bbm_pca100;
Rsq2_pca100 = 1 - sum((NETd - ybmpca100).^2)/sum((NETd).^2)

bbm_pca1571=bmpca1571\NETd;
ybmpca1571=bmpca1571*bbm_pca1571;
Rsq2_pca1571 = 1 - sum((NETd - ybmpca1571).^2)/sum((NETd).^2)
