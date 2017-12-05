# NOTES & Things to discuss in the meeting


## 05/12
BM visualisation:
1. For group Znet (full correlation), I get 'inf' for loads of values. 
 * In *nets_groupmean*, the way of converting t-score to z-score needs explaination!
 * I replace all 'inf' with the largest non-inf values, '-inf' with the smallest non-inf value. 



## 12/11
1. How do we treat outliers?
* Exclude variables with outliers
* Make the outliers into 0 but keep the variable.

2. If we exclude all variables with outliers, total number of SM variable gone down to 197 from 213. SDR SM gone down to 56 from 60 with alcohol 12 to 10, cognition 14 to 13, sensory 3 to 2 


## Week b/w 30/10
1. after sign-flipping some varaible like THC originally had strong correltion, the correlation gone weak

## Week b/w 13th Oct
1. SDR on unflipped data gives 60dim. 
* With sub-dims: 10	8	1	2	14	3	3	1	7	4	3	1	3
* Compared with 58dim. flipped dataset: 10 8 1 2 14 3 3 1 7 2 3 1 3
* 2 dims difference in physical health....

2. SDR on all unflipped 213 SDR variables gives 51 dim. as well.  (see error curve down the page)


## Week b/w 02/10/2017
1. In Steve's work, for SM, he first of all removed all ill-conditioned variables, then de-confound and calculate left eigenvectors.

2. For my SDR work:
* First I tried with de-confound within each domain after principal components been calculated. (SDR -> Badvars -> PC -> Deconf)
  * SDR SM dim. 47
* Then I switched to de-confound everything first and then go through SDR. (Deconf -> SDR ...)
  * SDR SM dim. 70
  * This cuased some ill-conditioned variables (e.g. with too many missing values or have too many same values) been taken into analysis.
  (Is this a good way of doing things?)
* To strictly follow Steve's order of method, I decided to de-confound just after removing ill-conditioned variables within each SDR domain.
(SDR -> Badvars -> Deconf -> PC)
  * SDR SM dim. 58
  
3. Factor rotation:
* We rotate structral correlation (canonical loadings), instead of rotating canoincal variables and then correlate with the observable variables.
  * rotate canoincal weights are not recommended

### Results
1. Visualisation:
* Showing 3d graph, use 'PCA_CCA' in workspace

2. 2-layer CV on all SDR variables 213: gives 51 dim.
* Error curve:  
![alt text](https://github.com/lzdh/SDR-CCA/blob/master/error_allSDRvars.jpg)



### Questions
1. How to do factor rotation in CCA??
* Do we rotate weights in SM and BM seperately? (using different tranformation matrices)
* if we are only interested in explaining the SM weights, can we just rotate the SM weights and compare with the unrotated BM weights? Or we rotate BM with the same transformation matrix used in SM?

2. It doesnt make sense to rotate canonical variables?

3. Buster stopped working... after I tried to add public key pairs with go password-less and it also failed...

4. Access to UK Biobank
