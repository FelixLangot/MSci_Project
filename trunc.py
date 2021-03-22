#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:40:35 2020

@author: félix langot
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 09:33:13 2020

@author: félix langot
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
    return (DT)**2 * (1 / (SX * (10000/3600))) * (9.1093856E-31 * 299792458**2 / (T_e * 1000 * 1.60217662E-19))**2 * ((L_ee/1000000)/(4*np.pi*4*(2.72548*6.6524587158E-29)**2 * (1+z)**4)) * 1/(2*np.sqrt(np.pi)) * (Xi(beta)/r_c) / 3.09E+25

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

# Thresold for the normal law:  95% => 1.96 | 90% =>1.64 | 80% => 1.28 | 68% => 1
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
# Plot figures of the distributions
# 0 = No
# 1 = yes
figplot = 1

# Generate 4 draws N(0,1) of NumSim sample
vec1       = np.zeros(NumObs)
vec2       = np.ones(NumObs)
N1_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N2_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N3_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N4_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N5_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N6_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N7_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)

#==============================================================================
# DATA
#==============================================================================
# Redshift and its unc
z = np.array([0.35,0.348114,0.3582,0.4064,0.415,0.423513,0.4703,0.48,0.498492,0.55,0.561539,0.571,0.58,0.58105,0.597129,0.62,0.635236,0.67,0.721856,0.7234,0.73,0.792])
z_unc = 1e-50*np.ones(NumObs) #z / 15
# Delta T and its unc
DT = np.array([-0.00155430052805673,-0.000815181971670572,-0.000748105556051039,-0.000956458467503602,-0.00119665463025181,-0.00113058275001967,-0.000644220854989454,-0.0010040294904678,-0.000591420873210084,-0.000674615246955078,-0.000637364124669414,-0.000560812161950115,-0.000617247144838919,-0.0012225721460276,-0.000581858082967865,-0.000645152416303552,-0.000672584790280117,-0.000558837648523966,-0.000710091268711629,-0.000618522155212186,-0.000830919905182418,-0.00068203636675998])
DT_unc = np.array([9.59448457117068E-05,0.000112316148858629,8.57126626310302E-05,9.16816621102044E-05,0.00010186997931053,8.97373969862464E-05,0.000109319955917826,0.000123474447638424,8.99823894682166E-05,9.41165707705255E-05,0.000080849098207692,8.47187564609134E-05,0.000101939221473428,0.00018652841620266,9.83833269503857E-05,9.94758136226199E-05,0.000146495845155973,8.64788125049503E-05,0.000168787636487639,0.000118087183591956,8.93893874900851E-05,8.52975141548692E-05])
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

plt.errorbar(np.arange(22), DT, yerr=DT_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('$D_T$')
plt.show() 
plt.errorbar(np.arange(22), beta, yerr=beta_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('β')
plt.show() 
plt.errorbar(np.arange(22), SX, yerr=SX_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('$SX$')
plt.show() 
plt.errorbar(np.arange(22), r_c, yerr=r_c_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('$r_c$')
plt.show() 
plt.errorbar(np.arange(22), Te, yerr=Te_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('$Te$')
plt.show() 
plt.errorbar(np.arange(22), L_ee, yerr=L_ee_unc,fmt='d',color='darkblue')
plt.xlabel('cluster')
plt.ylabel('$\Lambda$')
plt.show() 

#==============================================================================
# Distributions of the data with their measurement errors
#==============================================================================    
# Initialization of the stochastic parameters
z_dist    = np.zeros((NumObs,NumSim))
Te_dist   = np.zeros((NumObs,NumSim))
L_ee_dist = np.zeros((NumObs,NumSim))
SX_dist   = np.zeros((NumObs,NumSim))
beta_dist = np.zeros((NumObs,NumSim))
r_c_dist  = np.zeros((NumObs,NumSim))
DT_dist   = np.zeros((NumObs,NumSim))

z_dist0    = np.zeros((NumObs,NumSim))
Te_dist0   = np.zeros((NumObs,NumSim))
L_ee_dist0 = np.zeros((NumObs,NumSim))
SX_dist0   = np.zeros((NumObs,NumSim))
beta_dist0 = np.zeros((NumObs,NumSim))
r_c_dist0  = np.zeros((NumObs,NumSim))
DT_dist0   = np.zeros((NumObs,NumSim))

for ii in np.arange(NumObs):
    # 1
    XX1             = get_truncated_normal(mean=z[ii],    sd=z_unc[ii],    low=0, upp=2)
    z_dist[ii,:]    = XX1.rvs(size=NumSim,random_state=1)
    # 2
    XX2             = get_truncated_normal(mean=Te[ii],   sd=Te_unc[ii],   low=0, upp=15)
    Te_dist[ii,:]   = XX2.rvs(size=NumSim,random_state=2)
    XX3             = get_truncated_normal(mean=L_ee[ii], sd=L_ee_unc[ii], low=0, upp=1e-14)
    L_ee_dist[ii,:] = XX3.rvs(size=NumSim,random_state=2)
    # 3
    XX4             = get_truncated_normal(mean=SX[ii],   sd=SX_unc[ii],   low=1, upp=1000)
    SX_dist[ii,:]   = XX4.rvs(size=NumSim,random_state=3)
    XX5             = get_truncated_normal(mean=beta[ii], sd=beta_unc[ii], low=0, upp=2)
    beta_dist[ii,:] = XX5.rvs(size=NumSim,random_state=3)
    XX6             = get_truncated_normal(mean=r_c[ii],  sd=r_c_unc[ii],  low=2, upp=70)
    r_c_dist[ii,:]  = XX6.rvs(size=NumSim,random_state=3)
    # 4
    XX7             = get_truncated_normal(mean=DT[ii],   sd=DT_unc[ii],   low=-1, upp=-1e-6)
    DT_dist[ii,:]   = XX7.rvs(size=NumSim,random_state=4)

# Generate 4 draws N(0,1) of NumSim sample
vec1       = np.zeros(NumObs)
vec2       = np.ones(NumObs)
N1_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N2_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N3_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N4_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)
N5_dist    = gen_dist(vec1,vec2*Tnorm, NumSim)

for ii in np.arange(NumObs):
    # Group z
    z_dist0[ii,:]    = z[ii]    + N1_dist[ii,:]*(z_unc[ii]*Tnorm)    
    # Group (Te, Lambda)
    Te_dist0[ii,:]   = Te[ii]   + N2_dist[ii,:]*(Te_unc[ii]*Tnorm)
    L_ee_dist0[ii,:] = L_ee[ii] + N2_dist[ii,:]*(L_ee_unc[ii]*Tnorm)
    # Group (SX, beta, r_c)
    SX_dist0[ii,:]   = SX[ii]   + N3_dist[ii,:]*(SX_unc[ii]*Tnorm)
    beta_dist0[ii,:] = beta[ii] + N3_dist[ii,:]*(beta_unc[ii]*Tnorm)
    r_c_dist0[ii,:]  = r_c[ii]  + N3_dist[ii,:]*(r_c_unc[ii]*Tnorm)    
    # Group DT
    DT_dist0[ii,:]   = DT[ii]   + N4_dist[ii,:]*(DT_unc[ii]*Tnorm)


if figplot == 1:
    plt.hist(z_dist[1,:],100)
    plt.xlabel('z')
    plt.ylabel('Number of realization')
    plt.show()
    plt.hist(Te_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(Te_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('Te')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()
    plt.hist(L_ee_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(L_ee_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('L_ee')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()
    plt.hist(SX_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(SX_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('SX')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()
    plt.hist(beta_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(beta_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('beta')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()
    plt.hist(r_c_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(r_c_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('r_c')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()
    plt.hist(DT_dist[1,:],100, alpha=0.5, label='Trunc')
    plt.hist(DT_dist0[1,:],100, alpha=0.5, label='Raw')
    plt.xlabel('DT')
    plt.ylabel('Number of realization')
    plt.legend(loc='upper right')
    plt.show()

# Simulate the sample
D_a_dist   = D_a_eq(DT_dist, SX_dist, Te_dist, z_dist, beta_dist, r_c_dist, L_ee_dist)
D_a_dist_c = D_a_eq(DT_dist, SX_dist, Te_dist, z_dist, beta_dist, r_c_dist, L_ee_dist)

D_a_dist0   = D_a_eq(DT_dist0, SX_dist0, Te_dist0, z_dist0, beta_dist0, r_c_dist0, L_ee_dist0)
D_a_dist0_c = D_a_eq(DT_dist0, SX_dist0, Te_dist0, z_dist0, beta_dist0, r_c_dist0, L_ee_dist0)


# Prints each mean/error
errDalist_r   = []
meanDalist_r  = []
errDalist0_r  = []
meanDalist0_r = []

print("==========================================================")
print('Statistics for the raw simulated data')
for i in range(len(z)):
    D_a_mean_r   = np.mean(D_a_dist[i,:])
    D_a_error_r  = np.std(D_a_dist[i,:])
    D_a_mean0_r  = np.mean(D_a_dist0[i,:])
    D_a_error0_r = np.std(D_a_dist0[i,:])
    print('Mean D_a Truncatured O      = ', D_a_mean_r)
    print("Error on D_a  Truncatured O = ", D_a_error_r)
    print("----------------------------------------------------------")
    print('Mean D_a raw O              = ', D_a_mean0_r)
    print("Error on D_a  raw O         = ", D_a_error0_r)    
    print("==========================================================")
    errDalist_r.append(D_a_error_r)
    meanDalist_r.append(D_a_mean_r)
    errDalist_r.append(D_a_error0_r)
    meanDalist_r.append(D_a_mean0_r)

data_name = ['Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9','Cluster 10','Cluster 11',
             'Cluster 12','Cluster 13','Cluster 14','Cluster 15','Cluster 16','Cluster 17','Cluster 18','Cluster 19','Cluster 20','Cluster 21','Cluster 22',]
# To exclude the extreme values for each simulated data
LimH      = np.zeros(NumObs)
LimL      = np.zeros(NumObs)
LimH0     = np.zeros(NumObs)
LimL0     = np.zeros(NumObs)

for i in np.arange(NumObs):
    exclude0  = []
    exclude00 = []
    # exclude all values lower than 0
    aa      = np.where(D_a_dist[i,:]<0)
    aal     = aa[0].tolist()
    if not aal:
        D_a_dist_c[i,:]
    if aal:
        exclude0 = aal
        D_a_dist_c[i,exclude0] = np.nan
    
    aa0     = np.where(D_a_dist0[i,:]<0) 
    aal0     = aa0[0].tolist()
    if not aal0:
        D_a_dist0_c[i,:]
    if aal0:
        exclude00 = aal0
        D_a_dist0_c[i,exclude00] = np.nan

    if figplot ==1:
        qq = D_a_dist_c[i,:]
        ww = D_a_dist0_c[i,:]
        
        plt.hist(D_a_dist[i,:],100, alpha=0.5, label='Raw (O trunc)')
        plt.hist(qq[~np.isnan(qq)],100, alpha=0.5, label='Trunc>0 (O trunc)')
        plt.xlabel('Data')
        plt.ylabel('Number of realization')
        plt.title(data_name[i])
        plt.legend(loc='upper right')
        plt.show()
        
        plt.hist(qq[~np.isnan(qq)],100, alpha=0.5, label='Trunc>0 (O trunc)')
        plt.hist(ww[~np.isnan(ww)],100, alpha=0.5, label='Trunc>0 (O notrunc)')
        plt.xlabel('Data')
        plt.ylabel('Number of realization')
        plt.title(data_name[i])
        plt.legend(loc='upper right')
        plt.show()

    GridHist     = NumSim-len(exclude0)
    rr           = D_a_dist_c[i,:]
    counti, bini = np.histogram(rr[~np.isnan(rr)],GridHist) 

    GridHist0      = NumSim-len(exclude00)
    rr0            = D_a_dist0_c[i,:]
    counti0, bini0 = np.histogram(rr0[~np.isnan(rr0)],GridHist0) 

    LimH[i] = NumSim-1
    LimL[i] = 0 
    ttt  = np.where(np.cumsum(counti)/GridHist<.025)
    tttl = ttt[0].tolist()
    if not tttl:
        tttl = 0
        print('mass point at the bottom')
    LimL[i] = bini[np.max(tttl)]
    sss  = np.where(np.cumsum(counti)/GridHist>.975)
    sssl = sss[0].tolist() 
    if not sssl:
        sssl = GridHist-1
        print('mass point at the top')
    LimH[i] = bini[np.min(sssl)]
    
    LimH0[i] = NumSim-1
    LimL0[i] = 0 
    ttt0     = np.where(np.cumsum(counti0)/GridHist0<.025)
    tttl0    = ttt0[0].tolist()
    if not tttl0:
        tttl0 = 0
        print('mass point at the bottom 0')
    LimL0[i] = bini[np.max(tttl0)]
    sss0  = np.where(np.cumsum(counti0)/GridHist0>.975)
    sssl0 = sss0[0].tolist() 
    if len(sssl0)==GridHist0:
        sssl0 = 1
    if not sssl0:
        sssl0 = GridHist0-1
        print('mass point at the top 0')
    LimH0[i] = bini[np.min(sssl0)]
        
    excludeH = np.where(D_a_dist_c[i,:]>LimH[i]) 
    D_a_dist_c[i,excludeH] = np.nan
    excludeH = np.where(D_a_dist_c[i,:]<LimL[i]) 
    D_a_dist_c[i,excludeH] = np.nan

    zzz  = np.where(D_a_dist0_c[i,:]>LimH0[i]) 
    zzzl = zzz[0].tolist()
    if zzzl:
        excludeH0 = zzzl
        D_a_dist0_c[i,excludeH0] = np.nan
    zzw  = np.where(D_a_dist0_c[i,:]<LimL0[i]) 
    zzwl = zzw[0].tolist()
    if zzwl:
        excludeH0 = zzwl
        D_a_dist0_c[i,excludeH0] = np.nan
    
    
    if figplot == 1:
        plt.hist(D_a_dist_c[i,:],100, alpha=0.5,  label='0.025<Trunc<0.95 (O trunc)')
        plt.hist(D_a_dist0_c[i,:],100, alpha=0.5, label='0.025<Trunc<0.95 (O notrunc)')
        plt.xlabel('Corrected data')   
        plt.ylabel('Number of realization')
        plt.title(data_name[i])
        plt.legend(loc='upper right')
        plt.show()

errDalist  = []
meanDalist = []   
errDalist0  = []
meanDalist0 = []    
print("==========================================================")
print('Statistics for the corrected simulated data')
for i in range(len(z)):
    D_a_mean  = np.nanmean(D_a_dist_c[i,:])
    D_a_error = np.nanstd(D_a_dist_c[i,:])
    D_a0_mean  = np.nanmean(D_a_dist0_c[i,:])
    D_a0_error = np.nanstd(D_a_dist0_c[i,:])
    print('Mean D_a Truncatured O      = ', D_a_mean)
    print("Error on D_a Truncatured O  =  ", D_a_error)
    print("----------------------------------------------------------")
    print('Mean D_a  Raw 0             = ', D_a0_mean)
    print("Error on D_a Raw 0          = ", D_a0_error)
    print("==========================================================")
    errDalist.append(D_a_error)
    meanDalist.append(D_a_mean)
    errDalist0.append(D_a0_error)
    meanDalist0.append(D_a0_mean)

