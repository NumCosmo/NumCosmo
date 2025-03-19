import numpy as np
from scipy.special import erf

def expected_value(mu, sd, c):

    cut = np.log(c)
    
    A = (cut - mu) / sd 
    
    B = ( 1 / np.sqrt(2 * np.pi) ) * np.exp( -0.5 * (A ** 2) )
    
    C = 1 - 0.5 * (  1 + erf( A / np.sqrt(2) )  )
    
    correction = (sd * B / C )
    
    return mu + correction


def standard_deviation(mu, sd, c):
    
    cut = np.log(c)
    
    A = (cut - mu) / sd 

    B = ( 1 / np.sqrt(2 * np.pi) ) * np.exp( -0.5 * ( A ** 2 ) )

    C = 1 - 0.5 * ( 1 + erf(A / np.sqrt(2)) )

    correction = np.sqrt( 1 + ( A * B / C  ) - ( B / C ) ** 2 )
        
    return sd * correction
