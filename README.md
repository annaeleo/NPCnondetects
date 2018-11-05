# NPCnondetects
The function Mv_twosample_perm_test.r performs a k-variate permutation test for two independent samples.
It requires the following input:

data : an array of (k + 1) columns. The first contain the groups and other the observations for each variable.
stat : specifies the test statistics: "dm" for difference of means; "ad" for Anderson-Darling; "ad.Ruym" for Anderson-Darling with Ruy
alt :  specifies if we want a two-sided test ("two.sided"), a one tailed test first_group > second_group ("greater") or a one tailed test first_group < second_group ("less)
f.comb : combination function. On default is "F" for "Fisher" combination function

@Output
it prints on screen:
- the p-value of the test;
- who is the first and the second sample;
- the alternative;
it returns a list of objects: 
- T.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the test statistics. Element in the first row is the observed test statistics. 
- P.distribution is a vector or array of (B+1)-elements i.e. the permutation distribution of the p-value statistics. Element in the first row is the observed p-value statistics.

#Example:
XY <- read.csv("Name of dataset.csv",header = TRUE, sep=";") #Replace "Name of dataset.csv" with your dataset name

source("Mv_twosample_perm_test.r")

#Mann-Whitney test
res.MW <- perm.data(data = XY, stat = "mw", B=1000, alt="two.sided")

#Anderson-Darling test
res.AD <- perm.data(data = XY, stat = "ad_Ruym", B=1000, alt="two.sided")

#Difference of means
res.DM <- perm.data(data = XY, stat = "dm", B=1000, alt="two.sided")

