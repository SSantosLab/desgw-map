import numpy as np
from scipy.interpolate import interp1d

#how area scale as number of event
def max_scale(N_max,alpha):
    return (2.*(1.-0.5**(1./N_max)))**(alpha/3.)

#solve for (p,N_max)    
def solve_p_N(area_left,rate,T_left, alpha, area_bar, area_bar_p):
    prange=np.arange(0,1,0.01)	#try a range of p
    N_max=np.zeros(np.size(prange))
    N_max_trial=np.arange(1.,50.,0.1)	#try a range of N_max

    # produce a useful version of the area_bar, area_bar_p curve
    area_bar=interp1d(area_bar_p,area_bar)	
    
    # solve the area_left function, 
    # find pairs of (p,N_max) that satisfy the function
    N_comp=max_scale(N_max_trial,alpha)*rate*T_left/N_max_trial	

    for i in range(0,prange.size) :
        pp = prange[i]
        ix=  np.argmin(abs(area_left-area_bar(pp)*N_comp))
        N_max[i]=N_max_trial[ix]
    
    #find the pair (p,N_max) that maximize the p_tot 
    arg=np.argmax(prange*rate*T_left/N_max)
    return prange[arg],N_max[arg]


#
# This is the average event- we evaluate the probabilty for the
# average event, so that we can consider the real event
# this is the decisison we would make for the average event
# (you would use this p on the average map to find the average area)
# 
#   inputs: 
#       area_left area that could be observed this season
#       rate      effective rate of triggers
#       T_left    time left in the season
#       cl, Ap    this is the cumulative  distribution function of area
#                   for a fixed probability p ( p_gw in our notation), 
#                   as measured on the sims
#       area_bar   this pair is the average area of the sims as a 
#       area_bar_p function of probability p (i.e. we vary p_gw, measure area_bar)
#       
def evaluate_average_event(area_left, T_left, rate, cl, Ap, area_bar, area_bar_p) :
    #
    # solve the general case
    #

    #the scaling index of the area, alpha=1.8 for 2015
    alpha = 1.8

    # solve the implicit equation for p and N_max
    p,N_max=solve_p_N( area_left, rate, T_left, alpha, area_bar, area_bar_p )	
    # N_max: how loud is the softest one 
    # p: the probability that maximizes total probability if all events were average

    return p, N_max

