## spinhelicity.sage
## This is my implementation of Pauli matrices, spinors, spin helicity 
##   using Sage Math
## The main reference 
##   Scattering Amplitudes
######################################################################################
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                  
##                                                     
## 20150417                                 
##                                                   
## This program, along with all its code, is free software; you can redistribute it 
## and/or modify it under the terms of the GNU General Public License as published by 
## the Free Software Foundation; either version 2 of the License, or 
## (at your option) any later version.            
##                                                   
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of        
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    
## GNU General Public License for more details.               
##                                                          
## You can have received a copy of the GNU General Public License   
## along with this program; if not, write to the Free Software Foundation, Inc.,  
## S1 Franklin Street, Fifth Floor, Boston, MA      
## 02110-1301, USA                                             
##                                                          
## Governing the ethics of using this program, I default to the Caltech Honor Code: 
## ``No member of the Caltech community shall take unfair advantage of    
## any other member of the Caltech community.''                       
##                                                      
## If you like what I'm doing and would like to help and contribute support,   
## please take a look at my crowdfunding campaign at ernestyalumni.tilt.com and Patreon 
## read my mission statement and give your financial support, no matter how small 
## or large, if you can and to keep checking my ernestyalumni.wordpress.com blog 
## and various social media channels for updates as I try to keep putting out 
## great stuff.                                                      
##                                                                
## Fund Science! & Help Ernest in Physics Research! : quantum super-A-polynomials 
## - researched by Ernest Yeung
##                                      
## ernestyalumni.tilt.com                                 
##                                                      
## Facebook     : ernestyalumni                      
## gmail        : ernestyalumni                      
## google       : ernestyalumni                               
## linkedin     : ernestyalumni                                       
## Patreon      : ernestyalumni
## Tilt/Open    : ernestyalumni                               
## tumblr       : ernestyalumni                                 
## twitter      : ernestyalumni                   
## youtube      : ernestyalumni                                       
## wordpress    : ernestyalumni              
##  
##                                     
######################################################################################  
## spinhelicity.py
## 
## We've always defined ourselves by the ability to overcome the impossible. 
## And we count these moments, these moments when we dared to aim higher, 
## to break barriers, to reach for the stars,to make the unknown known. 
## We count these moments as our proudest achievements, 
## but we've lost all that.
## And perhaps we've just forgotten that we are still pioneers, and we've barely 
## begun, and 
## our greatest accomplishments cannot be behind us, but our destiny lies above us. 
## -Interstellar (2014)
##
## My mission is to continue to advance collective learning and basic research in 
## theoretical and mathematical physics.  My efforts include opening up the 
## research process and research tools so that in the spirit of open-source software 
## and collaboration, others can benefit and help improve upon the process and tools.
##
## Facebook     : ernestyalumni  
## github       : ernestyalumni                                                 
## gmail        : ernestyalumni                                         
## linkedin     : ernestyalumni                                       
## tumblr       : ernestyalumni                                           
## twitter      : ernestyalumni                                     
## youtube      : ernestyalumni                                         
## 
## Ernest Yeung was supported by Mr. and Mrs. C.W. Yeung, Prof. Robert A. Rosenstone, 
## Michael Drown, Arvid Kingl, Mr. and Mrs. Valerie Cheng, and the 
## Foundation for Polish Sciences, Warsaw University.  
##
## This code is open-source, governed by the Creative Common license.  
## Use of code is governed by the Caltech Honor Code: 
## ``No member of the Caltech community shall take unfair advantage of 
## any other member of the Caltech community.'' 
##
######################################################################################

################################################################################
##### N : Set up 4-manifold N for 4-momentum k's
## see my Notes and Prospectives on Scattering Amplitudes 2015 on 
## ernestyalumni.wordpress.com for the rationale of the manifold N for 4-momentum
##

N = Manifold(4,'N',r'N',start_index=0)
U_k = N.open_subset('U_k')
cart_k.<k0,k1,k2,k3> = U_k.chart(r'k0:(-oo,oo):k^0 k1:(-oo,oo):k^1 k2:(-oo,oo):k^2 k3:(-oo,oo):k^3' )
sph_k.<k0sph,kr,kph,kth> = U_k.chart(r'k0sph:(-oo,oo):k^0 kr:(0,oo) kph:(0,2*pi) kth:(0,pi)')

ktransit_sph_to_cart = sph_k.transition_map( cart_k , [ k0sph, kr * sin(kth) *cos(kph) , kr * sin(kth) * sin(kph) , kr * cos(kth) ] )

ktransit_cart_to_sph = ktransit_sph_to_cart.set_inverse( k0 , sqrt( k1^2 + k2^2 + k3^2 ) , atan2(k2, k1) , atan2( sqrt( k1^2 + k2^2) , k3 ) )

#####
## Minkowski flat metric
#####

eta_k = N.lorentz_metric('eta_k')
eta_k[0,0], eta_k[1,1], eta_k[2,2], eta_k[3,3] = -1 , 1 , 1 , 1

eta_U_k = U_k.lorentz_metric('eta_U_k')
eta_U_k[0,0], eta_U_k[1,1], eta_U_k[2,2], eta_U_k[3,3] = -1 , 1, 1 , 1

##############################
########## END N setup
################################################################################


##############################
## test 4-momentum p,k
##############################

p0, p1, p2, p3 = var('p0 p1 p2 p3', domain="real")
assume(p0>0)
assume(p1,"real")
assume(p2,"real")
assume(p3,"real")
pr, pph, pth = var('pr pph pth', domain='real')
assume(pr>=0)
assume(pph,'real')
assume(pth,'real')

kk0, kk1, kk2, kk3 = var('kk0 kk1 kk2 kk3', domain='real')
assume(kk0>0)
assume(kk1,"real")
assume(kk2,"real")
assume(kk3,"real")
kkr, kkph, kkth = var('kkr kkph kkth', domain='real')
assume(kkr>=0)
assume(kkph,'real')
assume(kkth,'real')

ppt = N.point((p0,p1,p2,p3), cart_k, name='p_test')
kpt = N.point((kk0,kk1,kk2,kk3), cart_k, name='k_test')
pptsph = N.point((p0,pr,pph,pth), sph_k , name='p_sph_test')
kptsph = N.point((kk0,kkr,kkph,kkth), sph_k , name='k_sph_test')

# Examples of usage:
ppt.coord()
ppt.coord(sph_k)

##############################
########## END test p,k
##############################

################################################################################
## Pauli matrices

from sympy import LeviCivita
from sympy.physics.matrices import mgamma, msigma

def Commutation( A , B ):
    return ( A*B - B*A )

def AntiCommutation( A , B ):
    return ( A * B + B * A )

Pauli_m = [ matrix.identity(2) , matrix( 2 , [ 0 , 1 , 1 , 0 ] ) , matrix( 2 , [ 0 , -I , I , 0 ] ) , matrix( 2 , [ 1 , 0 , 0 , -1 ] ) ]

# gamma0 = Matrix( [ [ 0 , 0 , 1 , 0 ] , 
#                   [ 0 , 0 , 0 , 1 ] ,
#                   [ 1 , 0 , 0 , 0 ] ,
#                   [ 0 , 1 , 0 , 0 ] ] )  

gamma0 = block_matrix( [ [ 0 , matrix.identity(2) ] , [ matrix.identity(2) , 0 ] ] )
# follows Srednicki's convention, see (36.39), (38.25) # gamma0 is also beta matrix (38.8)

gammamu = [ gamma0 , ] + [ block_matrix( [ [ 0 , Pauli_m[i] ] , [ - Pauli_m[i] , 0 ] ] ) for i in range(1,4) ]

gamma_5 = I * prod( gammamu )
print "gamma_5: ", gamma_5

################################################################################

def pslash( p ):
    """
    pslash
    p = ( p0 , p1 , p2 , p3 )

    e.g. pslash( ppt.coord()[:] ) # or 
    pslash( ppt.coord() )
    """
    return sum( [ sum( [ eta_U_k[ nu , mu ].expr() * p[nu] * gammamu[mu] for mu in range(4) ] ) for nu in range(4) ] )


def phi_a( k ): # make phi_a twistor
    """
    phi_a - phi_a twistor
    phi_a = phi_a( k )
    k = ( omega_k , r_k , phi_k , theta_k )

    e.g. phi_a( pptsph.coord(sph_k) )
    """
    omega_k , r_k , phi_k , theta_k = k
    assume( omega_k > 0 )
    phi_a = sqrt( 2 * omega_k ) * Matrix( [ [ - sin( theta_k / 2 ) * exp( -I * phi_k ) ] , 
    	      	      		  [   cos( theta_k / 2 ) ] ] )
    return phi_a
				 

def make_square_ket( k ): # make square ket 
    """
    make_square_ket - make square ket 
    make_square_ket = make_square_ket( k )
    k = ( omega_k , phi_k , theta_k )
    """
    return matrix( 4 , phi_a( k ).list() + [0,0] )			 
				 
def make_angle_bra( k ): # make angle bra 
    """
    make_angle_bra - make angle bra 
    make_angle_bra = make_angle_bra( k )
    k = ( omega_k , phi_k , theta_k )
    """
    phi_dota = phi_a( k ).conjugate_transpose()
    return matrix( 1 , [0,0] + phi_dota.list() )

def make_angle_ket( k ):
    """
    make_angle_ket - make angle ket
    make_angle_ket = make_angle_ket( k )
    k = ( omega_k , phi_k , theta_k ) 
    """
    phidota =  Pauli_m[1] * phi_a( k )
    return matrix( 4 , [0,0] + phidota.list() )

def make_square_bra( k ):
    """
    make_square_bra - make square bra 
    make_square_bra = make_square_bra( k )
    k = ( omega_k , phi_k , theta_k ) 
    """
    phia = phi_a( k ).conjugate_transpose()
    return matrix( 1 , phia.list() + [0,0] )


