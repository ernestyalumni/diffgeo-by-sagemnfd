## spinheli.sage
## This is my implementation of Pauli matrices, spinors, spin helicity 
##   using Sage Math
## The main reference 
##   Scattering Amplitudes
################################################################################################
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                                                      
##                                                                                                            
## 20150410                                 
##                                                                                                            
## This program, along with all its code, is free software; you can redistribute it and/or modify             
## it under the terms of the GNU General Public License as published by                                       
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
## read my mission statement and give your financial support, no matter how small or large, if you can        
## and to keep checking my ernestyalumni.wordpress.com blog and various social media channels                 
## for updates as I try to keep putting out great stuff.                                                      
##                                                                                                            
## Fund Science! & Help Ernest in Physics Research! : quantum super-A-polynomials - researched by Ernest Yeung
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
################################################################################################  
## spinhelicity.py
## 
## We've always defined ourselves by the ability to overcome the impossible. And we count these 
## moments, these moments when we dared to aim higher, to break barriers, to reach for the stars, 
## to make the unknown known. We count these moments as our proudest achievements, 
## but we've lost all that.
## And perhaps we've just forgotten that we are still pioneers, and we've barely begun, and 
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
################################################################################################

load("Min4d_base.sage")

def quatre_p_make( fourptest , chart=cart ): # quatre = four
    """
    quatre_p_make = quatre_p_make( fourptest )
    fourptest = ( p0 , p1 , p2 , p3 )
    chart=cart (Sage Manifold chart) Other examples are chart=sph or chart=cyl 
    	       or chart=cart_U (the restriction on a domain with spherical and cylindrical coordinates)

    """
    return sum( [ fourptest[i] * chart.frame()[i] for i in range(4) ] )

########################################
# sanity check
##

p0test , p1test , p2test , p3test = var("p0test p1test p2test p3test" , domain="real" )
fourptest = ( p0test , p1test , p2test , p3test ) # four = quatre

ph_k , th_k , ph_p , th_p = var("ph_k th_k ph_p th_p" , domain="real")
omega_k , omega_p = var("omega_k omega_p" , domain="positive")
k_test = ( omega_k , ph_k , th_k )
p_test = ( omega_p , ph_p , th_p )

cartptest = quatre_p_make( ( t , x , y , z ) , cart_U)
sphPtest  = quatre_p_make( ( p0test , omega_p , 0 , 0 ) , sph )
cartPtest = quatre_p_make( ( omega_p , omega_p * sin( th_p ) * cos( ph_p ) , omega_p * sin( th_p ) * sin( ph_p ) , omega_p * cos( th_p ) ) , cart_U )

cartPtest.down(g_U)( cartPtest ).display() # - \omega^2 + \omega^2

sphPtest.display(  cart_U.frame() , sph )
cartptest.display( cart_U.frame() , sph )
########################################

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
    """
    return sum( [ sum( [ g_U[ nu , mu ].expr() * p[nu] * gammamu[mu] for mu in range(4) ] ) for nu in range(4) ] )

pslash(fourptest)

################################################################################
## twistor
################################################################################


twistor_1 = U.scalar_field( -sin( th / 2) * exp( -I * ph ) , chart=sph ) 
twistor_2 = U.scalar_field(  cos( th / 2 ) , chart=sph )

twistor_a = Matrix( [ [ twistor_1.expr(sph) ] , 
	    	      [ twistor_2.expr(sph) ] ] ) # Srednicki (50.9)

twistor_a * twistor_a.conjugate_transpose()
( twistor_a * twistor_a.conjugate_transpose() ).nrows()
( twistor_a * twistor_a.conjugate_transpose() ).ncols()

( twistor_a * twistor_a.conjugate_transpose() )[0,0].expand_trig( half_angles=True )
( twistor_a * twistor_a.conjugate_transpose() )[1,0].expand_trig( half_angles=True )
( twistor_a * twistor_a.conjugate_transpose() )[0,1].expand_trig( half_angles=True )
( twistor_a * twistor_a.conjugate_transpose() )[1,1].expand_trig( half_angles=True )

pslash(fourptest).subs( p0test == 1 ).subs( p3test == cos(th) ).subs( p1test == sin(th) * cos(ph) ).subs( p2test == sin(th) * sin(ph) )


def phi_a( k ): # make phi_a twistor
    """
    phi_a - phi_a twistor
    phi_a = phi_a( k )
    k = ( omega_k , phi_k , theta_k )
    """
    omega_k , phi_k , theta_k = k
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


