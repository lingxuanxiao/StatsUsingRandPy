# BayesPowerAnalysisVignettes: Vignettes for Power Analysis in Bayesian Data Analysis  #
#                                                                                      #
# Author: Kazuha_Yae                                                                   #
# Ref: Kruschke, J. K. (2015)Doing Bayesian Data Analysis A Tutorial with R, JAGS, and #
#          Stan(2nd Ed.). Elsevier Inc.                                                #
#                                                                                      #
# Program Version: Alpha 0.0.1                                                         #
# R Version: ≥ 4.1                                                                     #
# Depends: 'rstan', 'bayestestR', 'BayesFactor'                                        #
# License: GPL 3.0                                                                     #
# Bug Report: sulfonamides@163.com                                                     #
#======================================================================================#
# 0. check packages & import packages
pkgs = c('rstan', 'bayestestR', 'BayesFactor')
for(i in 1:length(pkgs)){if(!(pkgs[i] %in% installed.packages()[, 'Package'])){install.packages(pkgs[i])}}
for(i in 1:length(pkgs)) library(pkgs[i])

# 1. bayestestR 与 BayesFactor 的说明
# bayestestR与BayesFactor
# https://easystats.github.io/bayestestR/
# https://richarddmorey.github.io/BayesFactor/



#10.3758/PBR.16.2.225