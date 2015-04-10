## Schwarzchild.sage
## This is my implementation of Schwarzschild spacetime from Eric Gourgoulhon, Michal Bejger (2014) 
##   using these authors' sagemanifolds 0.7 package for Sage Math
## The main reference that I'll liberally copy from is the pdf SM_Schwarzschild on the authors' webpage 
##   for sagemanifolds
################################################################################################              
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                                                      
##                                                                                                            
## 20150408                                 
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

## 20150408 EY :  I am implementing a solution to 
## Tutorial 13 Schwarzschild
##   WE Heraeus International Winter School on Gravity and Light
##


M = Manifold(4, 'M', r'M' ) ; print M # default start index is 0 which is what you want

#####
## Physical constants (i.e. parameters)
#####
G_N = var('G_N') ; assume(G_N >= 0)
M_0 = var('M_0') ; assume(M_0 >= 0)


###############
## Charts on M
###############

U_sph  = M.open_subset('U_sph'  ) ; print U_sph
spher.<t,r,ph,th>   = U_sph.chart(r't:(-oo,oo) r:(0,+oo) ph:(0,2*pi):\phi th:(0,pi):\theta')

# U_cart = U_sph.open_subset('U_cart' ) ; print U_cart
# cart.<tc,x,y,z>     = U_cart.chart()

# transit_sph_to_cart = spher.transition_map(cart , [ t, r*sin(th)*cos(ph) , r*sin(th)*sin(ph) , r*cos(th) ] , restrictions1 = x^2+y^2+z^2!= 0, restrictions2=r>0)
# transit_cart_to_sph = transit_sph_to_cart.set_inverse( tc , sqrt( x^2 + y^2 + z^2 ) , atan2(y,x) , atan2( sqrt(x^2 + y^2), z )  )


###############
## Metric Structures
###############

# eta = U_cart.riemann_metric('eta', latex_name=r'\eta')
# eta[0,0] = -1
# eta[1,1] = 1
# eta[2,2] = 1
# eta[3,3] = 1

# M.set_default_chart(spher)
espher = spher.frame() # frame

# g = U_sph.riemann_metric('g',latex_name=r'g')
# g[espher,0,0] = - ( 1 - 2 * G_N * M_0 / r )
# g[espher,1,1] =  1 / ( 1 - 2 * G_N * M_0 / r )
# g[espher,2,2] = r**2 
# g[espher,3,3] = r**2 * (sin(th))**2

g = M.lorentz_metric('g')
g[0,0] = - ( 1 - 2 * G_N * M_0 / r )
g[1,1] =  1 / ( 1 - 2 * G_N * M_0 / r )
g[2,2] = r**2 
g[3,3] = r**2 * (sin(th))**2


#####
## Levi-Civita connection
#####
nab = g.connection(); print nab

#####
## Curvature
#####

R = g.riemann()

#####
## Ricci tensor
#####

Ric = g.ricci()

#####
## Curves in $M$
#####

tau = var("tau")

gamma0 = function('gamma0',tau)
gamma1 = function('gamma1',tau)
gamma2 = function('gamma2',tau)
gamma3 = function('gamma3',tau)
curve  = [ gamma0.diff(tau) , gamma1.diff(tau) , gamma2.diff(tau) , gamma3.diff(tau) ]

dotgamma = sum( [ curve[i] * espher[i] for i in range(4) ] ) 

g(dotgamma, dotgamma)


tp, rp, php, thp = var('tp rp php thp') # t' tprime, r' rprime, phi' phi prime, theta' theta prime

gammap = tp * espher[0] + rp * espher[1] + php * espher[2] + thp * espher[3]

L = g(gammap, gammap)

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Exercise 1 Question 1
##

# get out the coefficients you want clearly

L.expr().coefficients(tp)[1][0].factor().full_simplify()
L.expr().coefficients(rp)[1][0].factor().full_simplify()
L.expr().coefficients(php)[1][0].factor().full_simplify()
L.expr().coefficients(thp)[1][0].factor().full_simplify()

# END Tutorial 13 Exercise 1 Question 1


g(gammap, gammap ).expr().diff( tp )

g0p = gamma0.diff(tau)

L.expr().diff( tp ).full_simplify().subs( tp== g0p ).subs( r == gamma1 ).diff(tau).full_simplify()

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Exercise 1 Question 2
##

L.expr().diff(t) # 0
L.expr().diff(tp).full_simplify().subs( tp == g0p ).subs( r == gamma1).diff(tau).expand() # to take partial derivatives of partial derivatives, the partial derivative has to be a variable

nab.coef()[:]
print nab.coef()[0,1,0].expr().factor()

# END Tutorial 13 Exercise 1 Question 2

# Why stop here? Get the 3 other Euler-Lagrange equation applications giving us equations of motion.

#
# r 
# 
L.expr().diff( r )
L.expr().diff( rp ).full_simplify().subs( rp == gamma1.diff( tau ) ).subs( r == gamma1 ).diff( tau ).expand()

#
# \phi
# 
L.expr().diff( ph ) # 0 
L.expr().diff( php ).factor().subs( php == gamma2.diff( tau ) ).subs( r == gamma1 ).diff(tau)

#
# \theta
#

L.expr().diff( th ).factor() 
L.expr().diff( thp ).factor().subs( r == gamma1 ).subs( thp == gamma3.diff( tau ) ).subs( th == gamma3 ).diff(tau).factor()

	       
##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Schwarzschild Spacetime Exercise 1 Question 3
##

K_t = espher[0]
g.lie_der(K_t).display() # 0, as desired

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Schwarzschild Spacetime Exercise 1 Question 4
##

# the given conserved quantity $(K_t)_a(x^a)'(\lambda) = const.$ 
K_t.down(g)( dotgamma ).expr() # (2*G_N*M_0*D[0](gamma0)(tau) - r*D[0](gamma0)(tau))/r 

# I tried to plug directly into the equation of motion
# L.expr().diff(tp).full_simplify().subs( tp == g0p ).subs( r == gamma1 ).diff(tau).full_simplify().subs( gamma0.diff() == const ).subs( gamma0.diff().diff() == 0 ).subs( gamma1(tau) == r ).subs( gamma1(tau).diff() == rp ) # -4*G_N*M_0*const*rp/r^2

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Schwarzschild Spacetime Exercise 1 Question 5
##

X_1 = sin(ph) * espher[3] + cot(th) * cos(ph) * espher[2]
X_2 = cos(ph) * espher[3] - cot(th) * sin(ph) * espher[2]
X_3 = espher[2]

# EY : 20150409 What does it really mean a symmetry? Consider that the Lie derivative of g along X_3 ("z-direction") is 0
g.lie_der(X_3).display() # 0

# Lie derivative of g along X_1 or X_2, i.e. $\mathcal{L}_{X_1}g$, $\mathcal{L}_{X_2}g$ nontrivial

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Schwarzschild Spacetime Exercise 1 Question 5
##

X_3.down(g)( dotgamma).expr() # r^2*D[0](gamma2)(tau)
J = var('J')

##
# WE Heraeus Gravity and Light International Winter School 2015
# Tutorial 13 Schwarzschild Spacetime Exercise 1 Question 6
##

E = K_t.down(g)( dotgamma ).expr().subs( gamma0.diff() == tp ).collect(tp) # (2*G_N*M_0 - r)*tp/r
E_0 = var('E_0')
L_f = L.expr().subs( th = pi/2).collect( php ).subs( php == J/ r**2 )

L_f = L.expr().subs( th = pi/2 ).collect( php ).subs( php == J/ r**2 ).collect( tp ).subs( tp == sqrt(E_0) * r / ( 2*G_N*M_0 - r ) ).subs( thp == 0 ).full_simplify()

solve(L_f == -1 , E_0)[0].expand() # E_0 == rp^2 - 2*G_N*J^2*M_0/r^3 + 2*G_N*M_0/r + J^2/r^2 - 1  # L_f == -1 because I choose the East Coast convention for the metric.  # This is the answer


##############################
## Reissner-Nordstrom Geometry
##############################

Phi    = function('Phi', r )

Lambda = function('Lambda', r)

g_RN = M.lorentz_metric('g_RN')
g_RN[0,0] = - exp( 2 * Phi    )
g_RN[1,1] =   exp( 2 * Lambda )
g_RN[2,2] = r**2 
g_RN[3,3] = r**2 * (sin(th))**2
 
R_RN   = g_RN.riemann()

Ric_RN = g_RN.ricci()

# Symmetries:
g_RN.lie_der( K_t ).display() # 0 time directed time invariant
g_RN.lie_der( X_3 ).display() # 0 X_3 or "z" angular momentum is an invariant
