clear
clc

vars0=importdata('MF_psy_vars.mat');
vars=importdata('SDR_vars_PC.mat');
BMvars=importdata('SDR_NET1_PCnew.mat');
NET=importdata('NET.mat');
load('SMvars.mat')
SM_header=importdata('SM_header.txt');
header=importdata('headers.txt');


% Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
% EB=hcp2blocks('RESTRICTED_lzdh_3_20_2016_11_24_35.csv', [ ], false, vars0(:,1)); % change the filename to your version of the restricted file
% PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

% fam_info=readtable('Fam_info.txt');
% [C0,IA]=unique(fam_info.Family_ID);
% C1=hist(fam_info.Family_ID, unique(fam_info.Family_ID))';

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
SM_header=SM_header(varskeep');

% Normalise and deconfound the original SM
varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
% varsgrot(isnan(varsgrot))=0;

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu,ss,vv]=nets_svds(NETd,Nkeep); % SVD reduction
BMpc=uu*ss;

% Deconfound and demean SDR BM
% BMvars=nets_demean(BMvars-conf*(pinv(conf)*BMvars));   
% 
% %Apply PCA to further reduce the dimension of SDR BM
% [uu,ss,vv]=nets_svds(BMvars,100); % SVD reduction
% BMpc=uu*ss;

% Normalise and deconfound SDR SM
varsd=palm_inormal(vars);
for i=1:size(varsd,2)
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd(isnan(varsd))=0;

%%% CCA/
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(BMpc,varsd);

%%% CCA permutation testing
% grotRp=zeros(Nperm,Nkeep); clear grotRpval;
% for j=1:Nperm
%   j
%   [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(BMpc,varsd(PAPset(:,j),:));
% end
% for i=1:Nkeep;  % get FWE-corrected pvalues
%   grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
% end
% grotRpval
% Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%%% CCA weights for CCA mode 1
grotAAd = corr(grotU(:,1),NETd)'; % weights after deconfounding
grotBBd = corr(grotV(:,1),varsgrot,'rows','pairwise')'; % weights after deconfounding

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

% Alternatively
varsgrot(isnan(varsgrot))=0;
b1=grotV(:,1:3)\varsgrot;
yhat=grotV(:,1:3)*b1;
Rsq = 1 - sum((varsgrot - yhat).^2)/sum((varsgrot).^2)

b2=grotU(:,1:3)\NETd;
yhat2=grotU(:,1:3)*b2;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)


%%
header=readtable('Header_SDR.csv', 'ReadVariableNames',false);
header=table2cell(header);
figure(3)
subplot(3,1,1);
plot(grotB(:,1))
set(gca,'xtick',[1:55],'XTickLabel',[])
%set(gca,'xtick',[1,7,15,17,18,21,32,34,37,40,43,45,48,50,58])
subplot(3,1,2);
plot(grotB(:,2))
set(gca,'xtick',[1:55],'XTickLabel',[])
subplot(3,1,3);
plot(grotB(:,3))
set(gca,'xtick',[1:55],'XTickLabel',header,'XTickLabelRotation',90)

figure(4)
subplot(3,1,1);
plot(grotB(:,1),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',[],'XLim',[0,56])

%set(gca,'xtick',[1,7,15,17,18,21,32,34,37,40,43,45,48,50,58])
h2=subplot(3,1,2);
plot(grotB(:,2),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',[],'XLim',[0,56])

subplot(3,1,3);
plot(grotB(:,3),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',header(:,2),'XTickLabelRotation',90,'XLim',[0,56])


% Factor rotation
RW=rotatefactors(grotB(:,1:3),'Method','promax');

header=readtable('Header_SDR.csv', 'ReadVariableNames',false);
header=table2cell(header);
subplot(3,1,1);
plot(RW(:,1),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',[],'XLim',[0,55])

%set(gca,'xtick',[1,7,15,17,18,21,32,34,37,40,43,45,48,50,58])
h2=subplot(3,1,2);
plot(RW(:,2),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',[],'XLim',[0,55])

subplot(3,1,3);
plot(RW(:,3),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',header(:,2),'XTickLabelRotation',90,'XLim',[0,55])
