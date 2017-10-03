# Things to discuss in the meeting

## Week b/w 02/10/2017
1. In Steve's work, for SM, he first of all removed all ill-conditioned variables, then de-confound and calculate left eigenvectors.

2. For my SDR work:
* First I tried with de-confound within each domain after principal components been calculated. (SDR -> Badvars -> PC -> Deconf)
  * SDR SM dim. 55
* Then I switched to de-confound everything first and then go through SDR. (Deconf -> SDR ...)
  * SDR SM dim. 70
  * This cuased some ill-conditioned variables (e.g. with too many missing values or have too many same values) been taken into analysis.
  (Is this a good way of doing things?)
* To strictly follow Steve's order of method, I decided to de-confound just after removing ill-conditioned variables within each SDR domain.
(SDR -> Badvars -> Deconf -> PC)
  * SDR SM dim. 58
