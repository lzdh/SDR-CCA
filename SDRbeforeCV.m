clear
clc

vars0=importdata('MF_psy_vars.mat');
vars=importdata('SDR_vars_PC.mat');
BMvars=importdata('SDR_NET1_PC.mat');
NET=importdata('NET.mat');
load('SMvars.mat')
SM_header=importdata('SM_header.txt');
header=importdata('headers.txt');


% Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
% EB=hcp2blocks('RESTRICTED_lzdh_3_20_2016_11_24_35.csv', [ ], false, vars0(:,1)); % change the filename to your version of the restricted file

fam_info=readtable('Fam_info.txt');
[C0,IA]=unique(fam_info.Family_ID);
C1=hist(fam_info.Family_ID, unique(fam_info.Family_ID))';

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


% Normalise and deconfound the original SM
varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
varsgrot(isnan(varsgrot))=0;

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean

% Deconfound and demean SDR BM
BMvars=nets_demean(BMvars-conf*(pinv(conf)*BMvars));   

%Apply PCA to further reduce the dimension of SDR BM
[uu,ss,vv]=nets_svds(BMvars,Nkeep); % SVD reduction
BMpc=uu*ss;


% Normalise and deconfound SDR SM
varsd=palm_inormal(vars);
for i=1:size(varsd,2)
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd(isnan(varsd))=0;

Rsq_cvsm=[];
Rsq_cvbm=[];
%%%CV 
for k=1:50
    
    fprintf('\n Simulation # %d',k);
    
    % Split the whole cohort into 2 halves without breaking family structure
    C=[C0 C1];
    C=sortrows(C,-2);
    C(:,3)=zeros;
    for i=1:size(C,1)
       if rand(1)>0.2   %Change here to change the proportion of training set!
           C(i,3)=1;
       else C(i,3)=2;
       end
    end

    fam_info.group_ID=zeros(815,1);
    for i=1:size(fam_info,1)
        index=find(fam_info.Family_ID(i)==C(:,1));
        fam_info.group_ID(i)=C(index,3);
    end

    In=fam_info(fam_info.group_ID==1,1);
    Out=fam_info(fam_info.group_ID==2,1);
    In=table2array(In);
    Out=table2array(Out);
    
    % Generate held-in and held-out SM and BM sets
    index1=0;
    index2=0;
    for i=1:size(In,1)
        index1(i)=find(In(i)==vars0(:,1));
    end
    varsd1=varsd(index1',:); %the held-in set
    uu1=BMpc(index1',:);
    
    for i=1:size(Out,1)
        index2(i)=find(Out(i)==vars0(:,1));
    end
    varsd2=varsd(index2',:); %the held-out set
    uu1_test=BMpc(index2',:);
    
    %%% CCA
    [grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(uu1,varsd1); %CCA on held-in set
    [grotA2,grotB2,grotR2,grotU2,grotV2,grotstats2]=canoncorr(uu1_test,varsd2); %CCA on held-out set
    
    R_in(k,:)=grotR;
    R_out(k,:)=grotR2;
    
    R_A(k)=corr(grotA(:,1),grotA2(:,1));
    R_B(k)=corr(grotB(:,1),grotB2(:,1));
    
    %%% CV on CCA
    %Apply canonical loadings from training set to test set to get'psuedo-canonical variates'
    grotV_CV=varsd2*grotB; % psuedo-canonical variates for SM
    grotU_CV=uu1_test*grotA; % psuedo-canonical variates for NET
  
    for j=1:size(grotU_CV,2)
        R_CV(k,j)=corr(grotV_CV(:,j), grotU_CV(:,j));  %(1)
    end
    
    
    %%%Permutation test
    %prepare permutation matrix
%     EB_in=EB(index1,:);
%     EB_out=EB(index2,:);
%     PAPset_in=palm_quickperms([ ], EB_in, Nperm);      % the final matrix of permuations
%     PAPset_out=palm_quickperms([ ], EB_out, Nperm);     % the final matrix of permuations
% 
%     %%% CCA permutation testing
%     grotRp=[]; clear grotRpval grotRpval_out;
%     grotRp_out=[];
%     grotRp_cv=[];
%     for j=1:Nperm
%       [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,varsd1(PAPset_in(:,j),:));
%       [~,~,grotRp_out(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1_test,varsd2(PAPset_out(:,j),:));
%       
%       %Permutation for CV
%        grotVr_CV=varsd2(PAPset_out(:,j),:)*grotBr; %'psuedo-canonical variates' for SM
%        grotUr_CV=uu1_test*grotAr;
%        for q=1:size(grotUr_CV,2)
%             grotRp_cv(j,q)=corr(grotVr_CV(:,q), grotUr_CV(:,q));  %(1)
%        end  
%     end
%     
%     for i=1:25;  % get FWE-corrected pvalues
%       grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
%       grotRpval_out(i)=(1+sum(grotRp_out(2:end,1)>=grotR2(i)))/Nperm;
%       grotRpval_cv(k,i)=(1+sum(max(grotRp_cv(2:end,:)')>=R_CV(k,i)))/Nperm;  
%     end
%     Ncca_in(k)=sum(grotRpval<0.05)  % number of FWE-significant CCA components
%     Ncca_out(k)=sum(grotRpval_out<0.05)  % number of FWE-significant CCA components
%     Ncca_cv(k)=sum(grotRpval_cv(k,:)<0.05)  % number of FWE-significant CCA components
 
    %%%CV for amount of variance between explained by canonical variates in original variables
    for j=1:3
        TVE_SM(j,:)=corr(grotV(:,j),varsgrot(index1,:),'rows','pairwise').^2; %in & in
        TVE_Net(j,:)=corr(grotU(:,j),NETd(index1,1:size(NET,2))).^2; %NETd(index1,1:size(NET,2))
        
        Y1=[];
        Y2=[];
       for i=1:size(varsgrot,2)
            B1=regress(varsgrot(index2,i),grotV_CV(:,[1:(j-1) (j+1):3]));
            Y1(:,i)=varsgrot(index2,i)-grotV_CV(:,[1:(j-1) (j+1):3])*B1;
       end
       for i=1:size(NETd(index2,1:size(NET,2)),2)
           B2=regress(NETd(index2,i),grotU_CV(:,[1:(j-1) (j+1):3]));
           Y2(:,i)=NETd(index2,i)-grotU_CV(:,[1:(j-1) (j+1):3])*B2;
       end
       TVE_SM_CV(j,:)=corr(grotV_CV(:,j),Y1,'rows','pairwise').^2; %in & in
       TVE_Net_CV(j,:)=corr(grotU_CV(:,j),Y2).^2; %NETd(index1,1:size(NET,2))
        
        TVE_SM_out(j,:)=corr(grotV2(:,j),varsgrot(index2,:),'rows','pairwise').^2; %in & out
        TVE_Net_out(j,:)=corr(grotU2(:,j),NETd(index2,1:size(NET,2))).^2; %NETd(index1,1:size(NET,2))
        
        %%%Redundency index
%         RI_SM(j,:)=corr(grotU(:,j),varsgrot(index1,:),'rows','pairwise').^2;
%         RI_Net(j,:)=corr(grotV(:,j),NETd(index1,1:size(NET,2))).^2;
        
         CanVar=(grotU(:,j)+grotV(:,j))./2;
         RI_SM(j,:)=corr(CanVar,varsgrot(index1,:),'rows','pairwise').^2;
         RI_Net(j,:)=corr(CanVar,NETd(index1,1:size(NET,2))).^2;

    end
    
    
    %%% An alterative way of calculating R-squared value between CV-canonical
    %%% variates and original datasets.
    b1=[];
    b2=[];
    yhat=[];
    yhat2=[];
    
    b1=grotV_CV(:,1:3)\varsgrot(index2,:);
    yhat=grotV_CV(:,1:3)*b1;
    Rsq_cvsm(k) = 1 - sum((varsgrot(index2,:) - yhat).^2)/sum((varsgrot(index2,:)-mean(varsgrot(index2,:))).^2)

    b2=grotU_CV(:,1:3)\NETd(index2,:);
    yhat2=grotU_CV(:,1:3)*b2;
    Rsq_cvbm(k) = 1 - sum((NETd(index2,:) - yhat2).^2)/sum((NETd(index2,:)-mean(NETd(index2,:))).^2)
    
    [ObyCV_net(k),CVbyO_net(k)]=space_spanned(grotU2(:,1:3),grotU_CV(:,1:3));
    [ObyCV_sm(k),CVbyO_sm(k)]=space_spanned(grotV2(:,1:3),grotV_CV(:,1:3));
   
    
    AVE_SM(k,:)=mean(TVE_SM,2,'omitnan');
    AVE_Net(k,:)=mean(TVE_Net,2,'omitnan');
   
    AVE_SM_CV(k,:)=mean(TVE_SM_CV,2,'omitnan');
    AVE_Net_CV(k,:)=mean(TVE_Net_CV,2,'omitnan');
    AVE_SM_out(k,:)=mean(TVE_SM_out,2,'omitnan');
    AVE_Net_out(k,:)=mean(TVE_Net_out,2,'omitnan');
    AVE_RI_Net(k,:)=mean(RI_Net,2,'omitnan');
    AVE_RI_SM(k,:)=mean(RI_SM,2,'omitnan');
    
end
  