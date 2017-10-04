
clear
clc

header=importdata('headers.txt');
vars0=importdata('Data/MF_psy_vars.mat');
%vars_new=zeros(size(vars0,1),1);

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25 479 480 481 486]) vars0(:,[265 266]).^(1/3) ]); % Gaussianise.  479 480 481 race variables
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

varsd=palm_inormal(vars0); % Gaussianise
for i=1:size(varsd,2) % deconfound ignoring missing data
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
vars0=varsd;

Psychiatric = [59:102];
Alcohol = ([115:148]);
%Alcohol=[116:118 120:121 123:124 126:127 129:130 132:148];
Tobacco =([149:184]);
%Tobacco= [150:153 155:156 158:159 161:162 164:165 167:168 170:184];
Drug = ([107:114, 185:194]);
Cognition = ([195:240]);
Emotion = ([241:264]); 
%new!!!
FamHis = [41:58];
FemHeal = [32:40 487];
Demo = [5 8:13]; % SES + handedness
PhyHeal = [14:18 20:26];
Motor = [458:464];
Personality = [465:469];
Sensory = [103:106 470:478 482:485];


sub_area={Psychiatric, Alcohol, Tobacco, Drug, Cognition, Emotion, FamHis, ...
          FemHeal, Demo, PhyHeal, Motor, Personality, Sensory};
vars_new_in=[];
pcaloadings=[]; 
header_keep=[];
for d=1:length(sub_area)
        fprintf('\n Estimating dimension for sub-area %d', d);
        
        badvars=[]; %exclude bad variables in the sub-area
        vars_in=vars0(:,sub_area{d});
        header_d=header(sub_area{d});
        for i=1:size(vars_in,2)
          Y=vars_in(:,i); grotKEEP=~isnan(Y);  
          grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?

          %change the creterion when # of subjects changes
          if d==8
              if (sum(grotKEEP)>=200) & (std(Y(grotKEEP))>0)& (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95);
              % the 3rd thing above is:  is the size of the largest equal-values-group too large?
                i=i;
              else
                  [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot];
                  badvars=[badvars i];
              end
          else
              if (sum(grotKEEP)>=400) & (std(Y(grotKEEP))>0)& (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95);
                  % the 3rd thing above is:  is the size of the largest equal-values-group too large?
                i=i; % do nothing
              else
                [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot];
                badvars=[badvars i];
              end
          end
        end

        varskeep=setdiff([1:size(vars_in,2)], badvars);  
        varsd1=palm_inormal(vars_in(:,varskeep)); % Gaussianise
        header_keep=[header_keep; header_d(varskeep)];
        
        
        %%% Calculate dim. (minimising  PRESS) for each of the sub-areas in held-in set using LOOCV!
        error2=0;
        for n=1:size(varsd1,1)
            
            Train = varsd1([1:n-1 n+1:end],:);
            Test = varsd1(n,:);
            
            grot=Train; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
            varsdCOV = (grot'*grot) ./ (grotI'*grotI);
            varsdCOV=nearestSPD(varsdCOV);
            [vv1,dd1]=eigs(varsdCOV,size(Train,2));
            
            if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end
            Test(isnan(Test))=0; % impute missing data as zeros

            for j=1:min(size(vv1,2),20)
                for i=1:size(Test,2)
                    proj = Test(:,[1:i-1 i+1:end])*pinv(vv1([1:i-1 i+1:end],1:j))'*vv1(:,1:j)'; 
                    err2(i) = Test(i) - proj(i); %Pseudoinverse
                end
                error2(n,j) = sum(err2(:).^2);
            end
        end
        error2 = sum(error2); 
        N_pca(d)=find(error2==(min(error2)));
        
        %%% Calculate eigen-components
        varsdCOV=zeros(size(varsd1,1));
        for i=1:size(varsd1,1)
          for j=1:size(varsd1,1)
              grot=varsd1([i j],:); grot=grot(:,sum(isnan(grot))==0); varsdCOV(i,j)=grot(1,:)*grot(2,:)';
          end
        end
        varsdCOV2=nearestSPD(varsdCOV);
        [uu,dd]=eigs(varsdCOV2,N_pca(d));  % SVD (eigs actually)
        
        
        %%%%%%%
        varsd_t=varsd1;
        varsd_t(isnan(varsd_t))=0;
        [u,s,v]=svd(varsd_t);
        %%%%%%%
        
     
        grot=varsd1; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
        varsdCOV1 = (grot'*grot);
        varsdCOV1=nearestSPD(varsdCOV1); % minor adjustment: project onto the nearest valid covariance matrix
        [vv,dd]=eigs(varsdCOV1,N_pca(d));  % SVD (eigs actually)
        %[vv,~]=eigs(varsdCOV1,size(varsdCOV1,2));  % SVD (eigs actually)
                   
        %if dd(1,1)<dd(2,2) vv=fliplr(vv); end 
        pcaloadings{d}=vv;
        vars_new_in=[vars_new_in uu*sqrt(dd)];
       
end
    
%%
vars0=importdata('MF_psy_vars.mat');
vars=importdata('SDR_vars_PC.mat');
NET=importdata('SDR_NET1_PC.mat');
NET0=importdata('NET.mat');

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

badvars=[];
for i=1:size(vars0,2)
  Y=vars0(:,i); grotKEEP=~isnan(Y);  
  %change the creterion when # of subjects changes
  if (sum(grotKEEP)<800) 
    i=i; % do nothing
    badvars=[badvars i];
  end
end

varskeep=setdiff([1:size(vars0,2)], badvars);  
vars0=vars0(:,varskeep);


varsgrot=palm_inormal(vars0);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end

NET1=nets_demean(NET0);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean

vars(isnan(vars))=0;

% r=[];
% for j=1:size(vars,2)
%     r(j,:)=corr(vars(:,j),varsgrot,'rows','pairwise').^2;
% end
% rsum=nansum(r');

vars(:,54)=[];
rho=partialcorri(varsgrot,vars,'rows','pairwise');
rsquare=rho.^2;
varexp=sum(nanmean(rsquare));


rho1=partialcorri(NETd,NET);
rsquare1=rho1.^2;
varexp_bm=nansum(mean(rsquare1));

%-------PCA variance explained------
load('pca100_smdata.mat')
load('pca55_smdata.mat')
rho=corr(varsgrot,uu2,'rows','pairwise');
rsquare=rho.^2;
varexp=nansum(mean(rsquare));


%% R-squared value
load('pca100_smdata.mat')
load('pca55_smdata.mat')
vars0=importdata('Data/MF_psy_vars.mat');
vars=importdata('SDR_vars_PC.mat');
NET=importdata('SDR_NET1_PC.mat');
NET0=importdata('NET.mat');
load('Data/SMvars.mat')
SM_header=importdata('SM_header.txt');

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
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
%SM_header=SM_header(varskeep');

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

vars(isnan(vars))=0;
varsgrot(isnan(varsgrot))=0;

b1=X2\varsgrot;
yhat=X2*b1;
Rsq = 1 - sum((varsgrot - yhat).^2)/sum((varsgrot).^2)


b2=NET\NETd;
yhat2=NET*b2;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)

b_pca100=uu2\varsgrot;
ypca100=uu2*b_pca100;
Rsq_pca100 = 1 - sum((varsgrot - ypca100).^2)/sum((varsgrot).^2)


%%

Psychiatric = [59:102];
Alcohol = ([115:148]);
%Alcohol=[116:118 120:121 123:124 126:127 129:130 132:148];
Tobacco =([149:184]);
%Tobacco= [150:153 155:156 158:159 161:162 164:165 167:168 170:184];
Drug = ([107:114, 185:194]);
Cognition = ([195:240]);
Emotion = ([241:264]); 
%new!!!
FamHis = [41:58];
FemHeal = [32:40 487];
Demo = [5 8:13]; % SES + handedness
PhyHeal = [14:18 20:26];
Motor = [458:464];
Personality = [465:469];
Sensory = [103:106 470:478 482:485];


sub_area={Psychiatric, Alcohol, Tobacco, Drug, Cognition, Emotion, FamHis, ...
          FemHeal, Demo, PhyHeal, Motor, Personality, Sensory};

SMvarsi=[];
for i=1:13
    SMvarsi=[SMvarsi sub_area{i}];
end
SMvarsi=sort(SMvarsi);
SMvars=vars0(:,SMvarsi);


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


SM_header=header(SMvarsi(:,varskeep)');

%%
PC=[];
for d=1:length(sub_area)
        fprintf('\n Estimating dimension for sub-area %d', d);
        
        badvars=[]; %exclude bad variables in the sub-area
        vars_in=vars0(:,sub_area{d});
        header_d=header(sub_area{d});
        for i=1:size(vars_in,2)
          Y=vars_in(:,i); grotKEEP=~isnan(Y);  
          grot=(Y(grotKEEP)-median(Y(grotKEEP))).^2; grot=max(grot/mean(grot));  % do we have extreme outliers?

          %change the creterion when # of subjects changes
          if d==8
              if (sum(grotKEEP)>=200) & (std(Y(grotKEEP))>0)& (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95);
              % the 3rd thing above is:  is the size of the largest equal-values-group too large?
                i=i;
              else
                  [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot];
                  badvars=[badvars i];
              end
          else
              if (sum(grotKEEP)>=400) & (std(Y(grotKEEP))>0)& (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95);
                  % the 3rd thing above is:  is the size of the largest equal-values-group too large?
                i=i; % do nothing
              else
                [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot];
                badvars=[badvars i];
              end
          end
        end

        varskeep=setdiff([1:size(vars_in,2)], badvars);  
        varsd1=palm_inormal(vars_in(:,varskeep)); % Gaussianise
        header_keep=[header_keep; header_d(varskeep)];
    
%%% Calculate eigen-components
        varsdCOV=zeros(size(varsd1,1));
        for i=1:size(varsd1,1)
          for j=1:size(varsd1,1)
              grot=varsd1([i j],:); grot=grot(:,sum(isnan(grot))==0); varsdCOV(i,j)=grot(1,:)*grot(2,:)';
          end
        end
        varsdCOV2=nearestSPD(varsdCOV);
        [uu,dd]=eigs(varsdCOV2,N_pca(2,d));  % SVD (eigs actually)
        
        
        %%%%%%%
        varsd_t=varsd1;
        varsd_t(isnan(varsd_t))=0;
        [u,s,v]=svd(varsd_t);
        %%%%%%%
        
     
        grot=varsd1; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
        varsdCOV1 = (grot'*grot);
        varsdCOV1=nearestSPD(varsdCOV1); % minor adjustment: project onto the nearest valid covariance matrix
        [vv,dd]=eigs(varsdCOV1,N_pca(2,d));  % SVD (eigs actually)
        %[vv,~]=eigs(varsdCOV1,size(varsdCOV1,2));  % SVD (eigs actually)
                   
        %if dd(1,1)<dd(2,2) vv=fliplr(vv); end 
        pcaloadings{d}=vv;
        PC=[PC uu*sqrt(dd)];
       
end


%% SDR on the whole SM dataset

header=importdata('headers.txt');
vars0=importdata('Data/MF_psy_vars.mat');
%vars_new=zeros(size(vars0,1),1);

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25 479 480 481 486]) vars0(:,[265 266]).^(1/3) ]); % Gaussianise.  479 480 481 race variables
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

SMvars=SMvars(:,varskeep);

varsd=palm_inormal(SMvars); % Gaussianise
for i=1:size(varsd,2) % deconfound ignoring missing data
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
vars0=varsd;

pcaloadings=[]; 

varsd1=varsd;
           
        %%% Calculate dim. (minimising  PRESS) for each of the sub-areas in held-in set using LOOCV!
        error2=0;
        for n=1:size(varsd1,1)
            
            Train = varsd1([1:n-1 n+1:end],:);
            Test = varsd1(n,:);
            
            grot=Train; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
            varsdCOV = (grot'*grot) ./ (grotI'*grotI);
            varsdCOV=nearestSPD(varsdCOV);
            [vv1,dd1]=eigs(varsdCOV,size(Train,2));
            
            if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end
            Test(isnan(Test))=0; % impute missing data as zeros

            for j=1:min(size(vv1,2),200)
                for i=1:size(Test,2)
                    proj = Test(:,[1:i-1 i+1:end])*pinv(vv1([1:i-1 i+1:end],1:j))'*vv1(:,1:j)'; 
                    err2(i) = Test(i) - proj(i); %Pseudoinverse
                end
                error2(n,j) = sum(err2(:).^2);
            end
        end
        error2 = sum(error2); 
        N_pca =find(error2==(min(error2)));
        
        %%% Calculate eigen-components
        varsdCOV=zeros(size(varsd1,1));
        for i=1:size(varsd1,1)
          for j=1:size(varsd1,1)
              grot=varsd1([i j],:); grot=grot(:,sum(isnan(grot))==0); varsdCOV(i,j)=grot(1,:)*grot(2,:)';
          end
        end
        varsdCOV2=nearestSPD(varsdCOV);
        [uu,dd]=eigs(varsdCOV2,N_pca(d));  % SVD (eigs actually)
        
        
        %%%%%%%
        varsd_t=varsd1;
        varsd_t(isnan(varsd_t))=0;
        [u,s,v]=svd(varsd_t);
        %%%%%%%
        
     
        grot=varsd1; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
        varsdCOV1 = (grot'*grot);
        varsdCOV1=nearestSPD(varsdCOV1); % minor adjustment: project onto the nearest valid covariance matrix
        [vv,dd]=eigs(varsdCOV1,N_pca(d));  % SVD (eigs actually)
        %[vv,~]=eigs(varsdCOV1,size(varsdCOV1,2));  % SVD (eigs actually)
                   
        %if dd(1,1)<dd(2,2) vv=fliplr(vv); end 
        pcaloadings{d}=vv;
        vars_new_in=[vars_new_in uu*sqrt(dd)];
       
