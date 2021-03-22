#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:27:12 2020

@author: felixlangot
"""

# This is a program to simualte errors for MSci Project

import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from numpy.polynomial.polynomial import polyfit
from scipy.stats import chi2
import statsmodels.api as sm
import seaborn as sns
from scipy.stats import truncnorm
sns.set(color_codes=True)

# Compute a Normal distribution over a bounded interval
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm( (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd )

# Compute gamma functions
def Xi(beta):
    return (gamma(1.5*beta)/(gamma((1.5*beta)-0.5)))**2 * (gamma((3*beta)-0.5)/gamma(3*beta))

# Equation to determine errors from
def D_a_eq(DT, SX, T_e, z, beta, r_c, L_ee):
    return (DT)**2 * (1 / (SX * (10000/3600))) * (9.1093856E-31 * 299792458**2 / (T_e * 1000 * 1.60217662E-19))**2 * ((L_ee/1000000)/(4*np.pi*2**2*(2.72548*6.6524587158E-29)**2 * (1+z)**4)) * 1/(2*np.sqrt(np.pi)) * (Xi(beta)/r_c) / 3.09E+25

# Generates n columns of normally distribtued values for each value of the array
# With a control on the generator of ramdom variables <=> the draws are always the same
def gen_dist(array, array_unc, nn): 
    np.random.seed(0)
    array_dist = np.random.normal(array[0], array_unc[0], nn)
    for i in range(1, len(array)):
        np.random.seed(i)
        newrow = np.random.normal(array[i], array_unc[i], nn)
        array_dist = np.vstack([array_dist, newrow])
    return array_dist

# This function calculates the cosmological integration factor which depends on z.
def hiya(x):
    return ( ((1+x)**2 * (1 + (0.3*x))) - (x*(2+x)*0.7) )**(-0.5)

#print('########################################################################')
#print('############################## First part ##############################')
#print('########################################################################')

# Thresold for the normal law
# 95% => 1.96 | 90% =>1.64 | 80% => 1.28 | 68% => 1
Tnorm = 1  
# Number of parameters
NumPar = 7
# Number of simulations
NumSim =  10000# 100000  # 1000000 
# Number of observations
NumObs = 22
# Number of independant observations
NumFree = 4
# Choose the data: 
# 0 = initial values 
# 1 = values computed April 10th
NumData = 1
# Choice the error type
# 0 = all errors are independant
# 1 = one measurement error for each group of independant observations
# => 4 groups (Te, Lambda) (SX, beta, r_c), z and DT
Error_type = 1 
# Choose a bounded distribution for the parameters
# 0 = unbounded
# 1 = bounded
Trunc = 1
# Plot figures of the distributions
# 0 = No
# 1 = yes
figplot = 0

# Generate 4 draws N(0,1) of NumSim sample
vec1       = np.zeros(22)
vec2       = np.ones(22)
N1_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N2_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N3_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N4_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N5_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N6_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N7_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)

# Redshift and its unc
z = np.array([0.35,0.348114,0.3582,0.4064,0.415,0.423513,0.4703,0.48,0.498492,0.55,0.561539,0.571,0.58,0.58105,0.597129,0.62,0.635236,0.67,0.721856,0.7234,0.73,0.792])
z_unc = z / 15
# Delta T and its unc
DT = np.array([-0.00155430052805673,-0.000815181971670572,-0.000748105556051039,-0.000956458467503602,-0.00119665463025181,-0.00113058275001967,-0.000644220854989454,-0.0010040294904678,-0.000591420873210084,-0.000674615246955078,-0.000637364124669414,-0.000560812161950115,-0.000617247144838919,-0.0012225721460276,-0.000581858082967865,-0.000645152416303552,-0.000672584790280117,-0.000558837648523966,-0.000710091268711629,-0.000618522155212186,-0.000830919905182418,-0.00068203636675998])
DT_unc = np.array([9.59448457117068E-05,0.000112316148858629,8.57126626310302E-05,9.16816621102044E-05,0.00010186997931053,8.97373969862464E-05,0.000109319955917826,0.000123474447638424,8.99823894682166E-05,9.41165707705255E-05,0.000080849098207692,8.47187564609134E-05,0.000101939221473428,0.00018652841620266,9.83833269503857E-05,9.94758136226199E-05,0.000146495845155973,8.64788125049503E-05,0.000168787636487639,0.000118087183591956,8.93893874900851E-05,8.52975141548692E-05])

if NumData == 0:
    # Beta and its unc
    beta = np.array([0.757964,0.718976,0.813157,0.810782,0.608069,0.557648,1.24219,0.690375,0.644117,0.781599,0.641377,0.943045,0.800767,0.590183,0.529259,0.61237,0.670853,0.766558,0.710763,0.55382,0.773329,1.02021])
    beta_unc = np.array([0.0497747,0.0437258,0.103342,0.133337,0.0142615,0.0176717,0.449146,0.0530611,0.0475828,0.0652756,0.0356317,0.148997,0.121666,0.0121917,0.0323047,0.0385831,0.0460435,0.101498,0.101507,0.00653156,0.111659,0.242764])    
    # S_X and its unc
    SX = np.array([176.795195530726,53.9083064516129,33.1936961538461,53.1816843826656,761.323660714286,113.112850812408,9.64443282674772,32.2290820668693,15.0632716535433,9.56680622009569,6.92811927272727,4.25671854545454,8.8244701723376,248.107804370447,24.0627649769585,2.80953483742535,35.4755285540705,7.07238493150685,2.04445714285714,86.272869096582,6.76774918032787,2.25395221843003])
    SX_unc = np.array([15.8798212290503,4.86360483870968,4.02220384615385,3.98195911692559,62.4720535714286,47.0213707533235,1.35983927051672,4.59382262001627,1.73326535433071,0.829114234449761,0.549078109090909,0.554991812865497,1.12307591692444,19.1118418314256,4.3619930875576,0.261560783012608,4.36973511543135,1.01709739726027,0.329428085106383,4.82065579350447,0.843260163934426,0.254435221843003])
    # r_c and its unc
    r_c = np.array([34.1717124,28.627266,46.0875096,41.4140508,9.2523552,22.9045188,66.10266,30.6930264,27.3035892,34.0870884,26.0589768,51.560124,42.3725652,8.4457704,14.2696728,25.7122152,20.8047612,31.9341948,31.289724,4.7809854,39.2456592,57.804096])
    r_c_unc = np.array([3.84321372,6.07833,8.4177264,8.5691148,0.68069676,2.48688288,22.2053376,4.24377552,3.99028236,4.31270964,2.72045496,9.399168,9.1127256,0.59322408,2.87800812,3.2774826,2.98504764,6.717276,7.5943644,0.241941492,8.0932032,14.4964356]) 
    # T_e and its unc
    Te = np.array([7.78014,6.06551,5.02654,7.6613,6.23959,5.93566,7.90022,8.96459,10.1339,8.07639,7.66745,6.9855,7.99919,6.17429,6.41017,6.57803,7.88763,6.01785,6.10179,5.08212,7.02854,10.1561])
    Te_unc = np.array([0.951414,0.579698,1.05384,0.965753,0.606944,0.601316,1.57905,1.70625,2.89344,0.846069,0.92777,1.40405,1.18775,0.52882,0.960366,0.672422,0.998255,0.747407,0.906228,0.247599,0.975405,1.65734])
    # Lambda and its unc
    L_ee = np.array([2.91697013092296E-15,2.99975426152968E-15,2.77840813718233E-15,2.66433368095109E-15,2.82664776993781E-15,3.08294675753574E-15,2.63197894237175E-15,2.75685592416849E-15,2.79253570146022E-15,2.67088226077674E-15,2.79584861499756E-15,2.71623747114726E-15,2.60184950321253E-15,2.73536837848511E-15,2.75759905538023E-15,2.80308479631623E-15,2.57726979991298E-15,2.46815692705906E-15,2.67943635402969E-15,2.67954325777263E-15,2.48921177376781E-15,2.42535512027552E-15])
    L_ee_unc = np.array([1.22137900332226E-16,1.24463843217837E-16,1.84580842154112E-16,8.43378271108534E-17,1.11112885330488E-16,1.61304819328662E-16,1.79525359114539E-16,1.3930719697191E-16,2.17499938051192E-16,1.03976107164004E-16,1.25246296424165E-16,1.3801841701109E-16,1.9282516398961E-16,9.89977727698239E-17,2.14357029516505E-16,1.27701192185867E-16,1.32102358348748E-16,1.12448251494614E-16,1.43620215460641E-16,7.10395619555591E-17,2.29856495280741E-16,9.5128558038985E-17])

if NumData == 1:
    # Beta and its unc
    beta = np.array([0.757964,0.718976,0.639373,1.01192,0.64635,0.472049,1.24219,0.690375,0.694249,0.699439,0.641377,0.664896,0.800767,0.62068,0.53608,0.61237,0.817283,0.766558,0.464444,0.581962,0.773329,1.02021])
    beta_unc = np.array([0.0497747,0.0437258,0.0450767,0.133337,0.0142615,0.0176717,0.449146,0.0530611,0.0548642,0.0779211,0.0356317,0.0467357,0.121666,0.0140217,0.0323047,0.0385831,0.0782782,0.101498,0.0448218,0.0079473,0.111659,0.242764])
    # S_X and its unc
    SX = np.array([176.7951955,53.90830645,36.81623077,46.08412101,719.4334821,113.1128508,9.644432827,32.22908207,14.73475394,10.62824641,6.928119273,4.001986909,8.824470172,235.3130073,24.68073733,2.809534837,32.78580802,7.072384932,2.854431611,80.64592256,6.76774918,2.253952218])
    SX_unc = np.array([15.87982123,4.863604839,4.463930769,3.981959117,62.47205357,47.02137075,1.359839271,4.59382262,1.684325197,1.110251196,0.549078109,0.485839766,1.123075917,16.92756712,4.361993088,0.261560783,4.369735115,1.017097397,0.575343465,4.820655794,0.843260164,0.254435222])
    # r_c and its unc
    r_c = np.array([34.1717124,28.627266,34.1302368,56.214936,10.4036352,11.2809204,66.10266,30.6930264,30.0408804,27.5993796,26.0589768,28.5721128,42.3725652,9.402612,14.4856116,25.7122152,26.5983564,31.9341948,13.8293328,5.4130332,39.2456592,57.804096])
    r_c_unc = np.array([3.84321372,6.07833,5.14755,8.5691148,1.3613886,2.48688288,22.2053376,4.24377552,4.30699752,4.90001004,2.72045496,3.59630844,9.1127256,0.61617588,2.87800812,3.2774826,4.050759,6.717276,4.24356396,0.25610568,8.0932032,14.4964356])
    # T_e and its unc    
    Te = np.array([7.78014,6.06551,5.02654,7.6613,5.43676,5.93566,7.90022,8.96459,10.1339,8.07639,7.66745,6.9855,7.99919,5.83822,5.101656667,6.57803,5.43698,6.01785,6.10179,5.08212,7.02854,10.1561])
    Te_unc = np.array([0.951414,0.579698,1.05384,0.965753,0.708637,0.601316,1.57905,1.70625,2.89344,0.846069,0.92777,1.40405,1.18775,0.487224,0.758313667,0.672422,0.72233,0.747407,0.906228,0.247599,0.975405,1.65734])
    # Lambda and its unc
    L_ee = np.array([2.91697013092296E-15,2.99975426152968E-15,2.77840813718233E-15,2.66433368095109E-15,2.85558640908637E-15,3.08294675753574E-15,2.63197894237175E-15,2.75685592416849E-15,2.79253570146022E-15,2.67088226077674E-15,2.79584861499756E-15,2.71623747114726E-15,2.60184950321253E-15,2.8002910963161E-15,2.87023837095416E-15,2.80308479631623E-15,2.61796399514733E-15,2.46815692705906E-15,2.67943635402969E-15,2.67954325777263E-15,2.48921177376781E-15,2.42535512027552E-15])
    L_ee_unc = np.array([1.22137900332226E-16,1.24463843217837E-16,1.84580842154112E-16,8.43378271108534E-17,1.23450573045525E-16,1.61304819328662E-16,1.79525359114539E-16,1.3930719697191E-16,2.17499938051192E-16,1.03976107164004E-16,1.25246296424165E-16,1.3801841701109E-16,1.9282516398961E-16,1.04034054992689E-16,2.82535885632563E-16,1.27701192185867E-16,1.69203720877678E-16,1.12448251494614E-16,1.43620215460641E-16,7.10395619555591E-17,2.29856495280741E-16,9.5128558038985E-17])

# Distributions of the data with their measurement errors
    
# Initialization of the stochastic parameters
z_dist    = np.zeros((NumObs,NumSim))
Te_dist   = np.zeros((NumObs,NumSim))
L_ee_dist = np.zeros((NumObs,NumSim))
SX_dist   = np.zeros((NumObs,NumSim))
beta_dist = np.zeros((NumObs,NumSim))
r_c_dist  = np.zeros((NumObs,NumSim))
DT_dist   = np.zeros((NumObs,NumSim))

if Trunc == 0:
    # Choice of the error type
    if Error_type == 0:    
        for ii in np.arange(NumObs):
        # Group z
            z_dist[ii,:]    = z[ii]    + N1_dist[ii,:]*(z_unc[ii]*Tnorm)    
        # Group (Te, Lambda)
            Te_dist[ii,:]   = Te[ii]   + N2_dist[ii,:]*(Te_unc[ii]*Tnorm)
            L_ee_dist[ii,:] = L_ee[ii] + N3_dist[ii,:]*(L_ee_unc[ii]*Tnorm)
        # Group (SX, beta, r_c)
            SX_dist[ii,:]   = SX[ii]   + N4_dist[ii,:]*(SX_unc[ii]*Tnorm)
            beta_dist[ii,:] = beta[ii] + N5_dist[ii,:]*(beta_unc[ii]*Tnorm)
            r_c_dist[ii,:]  = r_c[ii]  + N6_dist[ii,:]*(r_c_unc[ii]*Tnorm)    
        # Group DT
            DT_dist[ii,:]   = DT[ii]   + N7_dist[ii,:]*(DT_unc[ii]*Tnorm)
    
    if Error_type == 1:    
        for ii in np.arange(NumObs):
        # Group z
            z_dist[ii,:]    = z[ii]    + N1_dist[ii,:]*(z_unc[ii]*Tnorm)    
        # Group (Te, Lambda)
            Te_dist[ii,:]   = Te[ii]   + N2_dist[ii,:]*(Te_unc[ii]*Tnorm)
            L_ee_dist[ii,:] = L_ee[ii] + N2_dist[ii,:]*(L_ee_unc[ii]*Tnorm)
        # Group (SX, beta, r_c)
            SX_dist[ii,:]   = SX[ii]   + N3_dist[ii,:]*(SX_unc[ii]*Tnorm)
            beta_dist[ii,:] = beta[ii] + N3_dist[ii,:]*(beta_unc[ii]*Tnorm)
            r_c_dist[ii,:]  = r_c[ii]  + N3_dist[ii,:]*(r_c_unc[ii]*Tnorm)    
        # Group DT
            DT_dist[ii,:]   = DT[ii]   + N4_dist[ii,:]*(DT_unc[ii]*Tnorm)
        
if Trunc ==1:
    for ii in np.arange(NumObs):
        XX1             = get_truncated_normal(mean=z[ii],    sd=z_unc[ii],    low=0, upp=2)
        z_dist[ii,:]    = XX1.rvs(NumSim)
        XX2             = get_truncated_normal(mean=Te[ii],   sd=Te_unc[ii],   low=0, upp=15)
        Te_dist[ii,:]   = XX2.rvs(NumSim)
        XX3             = get_truncated_normal(mean=L_ee[ii], sd=L_ee_unc[ii], low=0, upp=1e-14)
        L_ee_dist[ii,:] = XX3.rvs(NumSim)
        XX4             = get_truncated_normal(mean=SX[ii],   sd=SX_unc[ii],   low=1, upp=1000)
        SX_dist[ii,:]   = XX4.rvs(NumSim)
        XX5             = get_truncated_normal(mean=beta[ii], sd=beta_unc[ii], low=0, upp=2)
        beta_dist[ii,:] = XX5.rvs(NumSim)
        XX6             = get_truncated_normal(mean=r_c[ii],  sd=r_c_unc[ii],  low=2, upp=70)
        r_c_dist[ii,:]  = XX6.rvs(NumSim)
        XX7             = get_truncated_normal(mean=DT[ii],   sd=DT_unc[ii],   low=-1, upp=0)
        DT_dist[ii,:]   = XX7.rvs(NumSim)

if figplot == 1:
    plt.hist(z_dist[1,:],100)
    plt.xlabel('$z$')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(Te_dist[1,:],100)
    plt.xlabel('$T_e$')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(L_ee_dist[1,:],100)
    plt.xlabel('$L_{ee}$')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(SX_dist[1,:],100)
    plt.xlabel('$S_{X0}$')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(beta_dist[1,:],100)
    plt.xlabel('$ β $')
    plt.ylabel('Number of realization')
    plt.savefig('Graphs/betadist.eps', dpi=200)
    plt.show()
    plt.hist(r_c_dist[1,:],100)
    plt.xlabel('$r_c$')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(DT_dist[1,:],100)
    plt.xlabel('$\Delta T$')
    plt.ylabel('Number of realization')
    plt.show()

# Simulate the sample
D_a_dist   = D_a_eq(DT_dist, SX_dist, Te_dist, z_dist, beta_dist, r_c_dist, L_ee_dist)
D_a_dist_c = D_a_eq(DT_dist, SX_dist, Te_dist, z_dist, beta_dist, r_c_dist, L_ee_dist)

# Prints each mean/error
errDalist_r  = []
meanDalist_r = []
print("==========================================================")
print('Statistics for the raw simulated data')
for i in range(len(z)):
    D_a_mean_r = np.mean(D_a_dist[i,:])
    D_a_error_r = np.std(D_a_dist[i,:])
    print('Mean D_a = ', D_a_mean_r)
    print("Error on D_a = ", D_a_error_r)
    print("==========================================================")
    errDalist_r.append(D_a_error_r)
    meanDalist_r.append(D_a_mean_r)

# To exclude the extreme values for each simulated data
LimH = np.zeros(NumObs)
LimL = np.zeros(NumObs)
for i in np.arange(NumObs):
    # exclude all values lower than 0
    exclude0 = np.where(D_a_dist[i,:]<0) 
    D_a_dist_c[i,exclude0] = np.nan
    if figplot ==1:
        plt.hist(D_a_dist[i,:],100)
        plt.xlabel('Raw data')
        plt.ylabel('Number of realization')
        plt.show()
    GridHist = NumSim
    counti, bini = np.histogram(D_a_dist[i,:],GridHist) 
    LimH[i] = NumSim-1
    LimL[i] = 0 
    LimH[i] = bini[int(np.max(np.where(np.cumsum(counti)/GridHist<.975)))]#<.975)))]
    excludeH = np.where(D_a_dist[i,:]>LimH[i]) 
    D_a_dist_c[i,excludeH] = np.nan
    if figplot == 1:
        plt.hist(D_a_dist_c[i,:],100)
        plt.xlabel('Corrected data')   
        plt.ylabel('Number of realization')
        plt.show()
    
# Prints each mean/error for corrected data 
# WARNING : It exists NaN
errDalist  = []
meanDalist = []    
print("==========================================================")
print('Statistics for the corrected simulated data')
for i in range(len(z)):
    D_a_mean  = np.nanmean(D_a_dist_c[i,:])
    D_a_error = np.nanstd(D_a_dist_c[i,:])
    print('Mean D_a = ', D_a_mean)
    print("Error on D_a = ", D_a_error)
    print("==========================================================")
    errDalist.append(D_a_error)
    meanDalist.append(D_a_mean)

# Explained data
D_a      = D_a_eq(DT,SX,Te,z,beta,r_c,L_ee)  
errD_a   = np.asarray(errDalist)

#print('########################################################################')
#print('############################# Second part ##############################')
#print('########################################################################')

# Theoretical model 
c     = 2.99792458 * 10**5  # in km/s
conv  = 1000
Ndata = z.shape[0]

# Loop to calculate the integration factor for each z in the list
Ktup  = [] 
Klist = [] 
errK  = []
for i in np.arange(0, Ndata):
    Ktup.append(integrate.quad(lambda x: hiya(x), 0, z[i]))
    (value, error) = Ktup[i]
    Klist.append(value)
    errK.append(error)
K = np.asarray(Klist)

# Regression to obtain a estimation for H
# Assumption: Data = (1/H)*( (c*K)/(1+z) ) <=> H = ( (c*K)/(1+z) ) / Data 
# Statistical model: Y = H+epsilon  where Y = exp( log( (c*K)/(1+z) ) - log(Data) )
Xvec    = ( (c*K)/(1+z) )/conv
DLog    = Xvec/D_a
# Exogenous variable : a constant
ZZ = np.ones(22)

resultols1 = sm.OLS(DLog,ZZ).fit()

print('==============================================================================')
print('Number of simulated data               = ',NumSim)
print('Data set (0 for initial, 1 for last)   = ', NumData)
print('Error type (0 for independant, 1 else) = ', Error_type) 
print('==============================================================================')

print('==============================================================================')
print('Results using Ordinary least squares method (linear regression)')
print('==============================================================================')

print(resultols1.summary())

H0    = resultols1.params
stdH0 = resultols1.HC1_se

Dfit    = np.zeros(Ndata)
Dfit[:] = ( 1/H0 )*( (c*K)/(1+z) )/conv
χ       = np.sum( (Dfit[:]-D_a[:])**2/(errD_a[:]**2) )  
Pval    = chi2.sf(χ, NumObs-1)

print('==============================================================================')
print('χ statistic = ', χ)
print('Reduced χ   = ', χ/(NumObs-1) )
print('P-Value     = ', 1-Pval, '   P_χ(0.95% , 21 df) = ',11.591)
print('==============================================================================')


plt.errorbar(z, D_a, yerr=errD_a,fmt='d',color='darkblue')
limlow  = (1/(H0-Tnorm*stdH0)) * (c*K)/(1+z)/conv
limhigh = (1/(H0+Tnorm*stdH0)) * (c*K)/(1+z)/conv
plt.fill_between(z,limlow,limhigh,color='r',alpha=0.4)
plt.xlabel('z')
plt.ylabel('$D_a$ (Gpc)')
plt.savefig('Graphs/linearmethod.eps', dpi=200)
plt.show() 

print('==============================================================================')
print('Results using nonlinear method (Minimum of the distance)')
print('==============================================================================')
# Length for the H0 grid
N        = 100
χ_v      = np.zeros((N, 1))
χ_c      = np.zeros((N, NumSim))
numchi_c = np.zeros(NumSim)
solH0_c  = np.zeros(NumSim)
a        = 40*conv   # in km/s/Gpc 
b        = 110*conv
H0_vec   = np.arange(a, b, (b-a)/N)
D_a_dist_C = np.zeros((NumObs,NumSim))
Dfit_v   = np.zeros((N,Ndata))

# Find the data corresponging the the similation j
for j in np.arange(NumSim):
    D_a_dist_C[:,j] = D_a_dist[:,j]

# This loop calculates the χ^2 for different H0 values
for j in np.arange(0, N):
    H0 = H0_vec[j]
    Dfit_v[j,:]  = ( 1/H0 )*(c*K)/(1+z)
    χ_v[j] = np.sum( (Dfit_v[j,:]-D_a[:])**2/(errD_a[:]**2) )  

minval = np.min(χ_v)
minpos = np.argmin(χ_v)

# Polynomial approx for the complete grid of H
m0, n0, o0, p0, q0 = polyfit(H0_vec/conv,χ_v,4)  # a fourth degree polynomial looks to fit the data best
Poly_est      = q0*(H0_vec/conv)**4 + p0*(H0_vec/conv)**3 + o0*(H0_vec/conv)**2 + n0*(H0_vec/conv) + m0
# plot over all H in [Hmin, Hmax]
plt.plot(H0_vec/conv, Poly_est, '-', color='red', label='Polyfit: P(x)')
plt.scatter(H0_vec/conv,χ_v,color='darkblue',s=20, label='Empirical values')
plt.xlabel('$H_0~( km \cdot s^{-1} \cdot Mpc^{-1} )$')
plt.ylabel('$\chi^2$')
plt.legend()
plt.savefig('Graphs/chibroad.eps', dpi=200)
plt.show()

# Find the minimum and the χ^2 solution
# => To increas precision, we reduce the grid of H around its minimal value found in the first step
# new grid
NL            = 1000000
H0L_vec       = np.arange(H0_vec[minpos]*.8, H0_vec[minpos]*1.2, H0_vec[minpos]*(1.2-0.8)/NL)
minposχ_v     = np.max(np.where(H0_vec<H0_vec[minpos]*.8))
maxposχ_v     = np.min(np.where(H0_vec>H0_vec[minpos]*1.2))
H0S_vec       = H0_vec[minposχ_v:maxposχ_v] #np.arange(H0_vec[np.argmin(χ_v)]*.8, H0_vec[np.argmin(χ_v)]*1.2, H0_vec[np.argmin(χ_v)]*(1.2-0.8)/NL)
# Polynomial approx for a shorter grid of H
m, n, o, p, q = polyfit(H0S_vec/conv,χ_v[minposχ_v:maxposχ_v],4)  # a fourth degree polynomial looks to fit the data best
Poly_estL     = q*(H0L_vec/conv)**4 + p*(H0L_vec/conv)**3 + o*(H0L_vec/conv)**2 + n*(H0L_vec/conv) + m
# Polynomial derivative
Poly_est_devL = 4*q*(H0L_vec/conv)**3 + 3*p*(H0L_vec/conv)**2 + 2*o*(H0L_vec/conv) + n
# Solution using derivative of the polynomial
solH0         = H0L_vec[np.argmin(np.abs(Poly_est_devL))] 
χ_min         = Poly_estL[np.argmin(np.abs(Poly_est_devL))]
χ_red         = χ_min/(NumObs-1)
Pval2         = chi2.sf(χ_min, NumObs-1)
# plot over a shorter range for H in [.8*Hhat, 1.2*Hhat]
plt.plot(H0L_vec/conv, Poly_estL, '-', color='red', label='Polyfit: P(x)')
plt.scatter(H0S_vec/conv,χ_v[minposχ_v:maxposχ_v],color='darkblue',s=20, label='Empirical values')
plt.xlabel('$H_0~( km \cdot s^{-1} \cdot Mpc^{-1} )$')
plt.ylabel('$\chi^2$')
plt.legend()
plt.savefig('Graphs/chinarrow.eps', dpi=200)
plt.show()
# plot of the derivative: to check the unicity of the zero
plt.plot(H0L_vec/conv, Poly_est_devL, '-', color='red', label='Derivative of Polyfit: dP(x)/dx')
plt.xlabel('$H_0~( km \cdot s^{-1} \cdot Mpc^{-1} )$')
plt.ylabel('$dP(x)/dx$')
plt.legend()
plt.savefig('Graphs/polyfitderivative.eps', dpi=200)
plt.show()

# Compute the distribution of the minimum value for H0
for ii in np.arange(NumSim):
    for j in np.arange(0, N):
        H0 = H0_vec[j]
        Dfit_v[j,:]  = ( 1/H0 )*(c*K)/(1+z)
        χ_c[j,ii] = np.sum( (Dfit_v[j,:]-D_a_dist_C[:,ii])**2/(errD_a[:]**2) )
    # Find the minimum
    numchi_c[ii] = np.argmin(χ_c[:,ii])
    solH0_c[ii]  = H0_vec[int(numchi_c[ii])]

H0_c_mean = np.mean(solH0_c)/conv
H0_c_med  = np.median(solH0_c)/conv
stdH0     = np.std(solH0_c/conv)

# Characteristic of the complete distribution of H0
counts, bins = np.histogram(solH0_c/conv,N)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('$H_0~( km \cdot s^{-1} \cdot Mpc^{-1} )$')
plt.ylabel('Number of realization')
plt.savefig('Graphs/distribution.eps', dpi=200)
plt.show()
# Find the 95% of H0 around its mean value
Lim       = np.max(np.where(np.cumsum(counts)/N<.96))
test1     = solH0_c/conv-bins[Lim]
itemindex = np.where(test1>0)
excludedH = solH0_c[itemindex]
minH0     = np.min(excludedH)/conv
maxH0     = np.max(excludedH)/conv

Conflow95 = solH0/conv-stdH0*1.96
Confhigh95 = solH0/conv+stdH0*1.96 

plusminus95 = solH0/conv - Conflow95

Conflow68 = solH0/conv-stdH0
Confhigh68 = solH0/conv+stdH0 

plusminus68 = solH0/conv - Conflow68

empconflow95  = solH0/conv - minH0
empconfhigh95 = maxH0 - solH0/conv

empconflow68  = solH0/conv - minH0*1.96
empconfhigh68 = maxH0/1.96 - solH0/conv

print('==============================================================================')
print('Estimated value for H (95%)   =', solH0/conv, '±',plusminus95)
print('Estimated value for H (68%)   =', solH0/conv, '±',plusminus68)
print('Standard deviation for H      =', stdH0)
print('Centered Conf. Interval (95%) =', Conflow95, Confhigh95)
print('Centered Conf. Interval (68%) =', Conflow68, Confhigh68)
print('Mean of H distribution        =', H0_c_mean)
print('Median of H distribution      =', H0_c_med )
print('Empirical Conf. Interval      =', minH0, maxH0)
print('Final H0 result (95%)         =', solH0/conv, '+',empconfhigh95, '-', empconflow95)
print('Final H0 result (68%)         =', solH0/conv, '+',empconfhigh68, '-', empconflow68)
print('χ                             =',χ_min  )
print('Reduced χ                     =',χ_red  )
print('P-Value                       =', 1-Pval2, ' P_χ(0.95% , 21 df) = ',11.591)
print('==============================================================================')

limlow2  = (1/(minH0)) * (c*K)/(1+z)/conv
limhigh2 = (1/(maxH0)) * (c*K)/(1+z)/conv

limlow3  = (1/(Conflow95)) * (c*K)/(1+z)/conv
limhigh3 = (1/(Confhigh95)) * (c*K)/(1+z)/conv

plt.errorbar(z, D_a, yerr=errD_a,fmt='d',color='darkblue')
plt.fill_between(z,limlow2,limhigh2,color='r',alpha=0.4)
plt.xlabel('z')
plt.ylabel('$D_a$ (Gpc)')
#plt.legend(loc='upper left') 
#plt.savefig('Graphs/nonlinear.eps', dpi=200)
plt.show() 

plt.errorbar(z, D_a, yerr=errD_a,fmt='d',color='darkblue')
plt.fill_between(z, limlow, limhigh, color='g', alpha=0.4, label='linear')
plt.fill_between(z, limlow2, limhigh2, color='r', alpha=0.4, label='nonlinearwrong')
plt.fill_between(z, limlow3, limhigh3, color='b', alpha=0.4, label='nonlinear95')
plt.xlabel('z')
plt.ylabel('$D_a$ (Gpc)')
plt.legend(loc='upper left') 
#plt.savefig('Graphs/Final.eps', dpi=200)
plt.show() 


plt.errorbar(z, D_a, yerr=errD_a, fmt='d', color='darkblue')
plt.fill_between(z, limlow2, limhigh2, color='r', alpha=0.4, label='nonlinear')
plt.fill_between(z, limlow, limhigh, color='g', alpha=0.4, label='linear')
plt.xlabel('z')
plt.ylabel('log($D_a$)')
plt.ticklabel_format(axis="y", style="plain", scilimits=(0,0))
plt.legend(loc='upper left')
plt.yscale("log")
#plt.savefig('Graphs/logFinal.eps', dpi=200)
plt.show()
