# BayesPowerAnalysis: Power Analysis in Bayesian Data Analysis                   #
#                                                                                #
# Author: Kazuha_Yae                                                             #
# Ref: Kruschke, J. K. (2015)Doing Bayesian Data Analysis A Tutorial with R, JAGS#
#          , and Stan(2nd Ed.). Elsevier Inc.                                    #
#                                                                                #
# Program Version: Alpha 0.0.1                                                   #
# R Version: ≥ 4.1                                                               #
# Depends: rstan, R6                                                             #
# License: GPL 3.0                                                               #
# Bug Report: sulfonamides@163.com                                               #
#================================================================================#
#-1. Check packages
if('rstan' %in% installed.packages()[, 'Package'])
if('R6'    %in% installed.packages()[, 'Package'])
