## S2.R
## This is my implementation of SM_sphere_S2 Sphere S^2 from Eric Gourgoulhon, Michal Bejger (2014) 
##   using these authors' sagemanifolds 0.7 package for Sage Math
## The main reference that I'll liberally copy from is the pdf SM_sphere_S2 on the authors' webpage 
##   for sagemanifolds
################################################################################################              
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                                                      
##                                                                                                            
## 20150322                                                                                                   
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

## 20150322 EY : This is v2, version 2, because I am implementing a solution to 
## Tutorial 11 Symmetry, Exercise 3: Lie derivative, 
##   WE Heraeus International Winter School on Gravity and Light
##

S2 = Manifold(2, 'S^2', r'\mathbb{S}^2', start_index=1) ; print S2

#
## Charts on S^2
#

# U := S^2\{N}
# V := S^2\{S}

U = S2.open_subset('U') ; print U
V = S2.open_subset('V') ; print V

S2.declare_union(U,V) # declare U, V cover S^2
U.union(V)

stereoN.<x,y> = U.chart()

S = U.point((0,0), name='S') ; print S

stereoS.<xp,yp> = V.chart("xp:x' yp:y'")

N = V.point((0,0), name='N') ; print N

trans = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)), intersection_name='W', restrictions1=x^2+y^2!=0, restrictions2=xp^2+yp^2!=0)

W = U.intersection(V)

inv_trans = trans.inverse()

stereoN_W = stereoN.restrict(W)
stereoS_W = stereoS.restrict(W)

#####
## Scalar fields
#####

f0 = S2.scalar_field({stereoN: function('F',x,y)}, name='f_0') ; f0.display()
f0.add_expr_by_continuation(stereoS, W) ; f0.display()

#####
## Vector fields
#####

u = U.vector_field('u')
u[:] = [function('u_x',x,y), function('u_y', x,y) ] ; u.display()

#####
## Embedding of S^2 into R^3
#####

R3 = Manifold(3, 'R^3', r'\mathbb{R}^3', start_index=1)
cart.<X,Y,Z> = R3.chart() ; cart

Phi = S2.diff_mapping(R3, {(stereoN, cart): [2*x/(1+x^2+y^2),2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)], (stereoS, cart): [2*xp/(1+xp^2+yp^2),2*yp/(1+xp^2+yp^2), (1-xp^2-yp^2)/(1+xp^2+yp^2)]}, name='Phi',latex_name=r'\Phi') ; Phi.display()

#####
## Riemannian metric on S^2
#####

h = R3.riemann_metric('h')
h[1,1], h[2,2], h[3,3] = 1,1,1 ; h.display()

g = S2.riemann_metric('g')
g.set( Phi.pullback(h) ) ; print g


#####
## Curvature
#####

Riem = g.riemann() ; print Riem
Riem.display()

Ric = g.ricci() ; Ric.display()

R = g.ricci_scalar() ; R.display()

#####
## Spherical coordinates
#####

A = W.open_subset('A', coord_def={stereoN_W: (y!=0, x<0)})
stereoN_A = stereoN_W.restrict(A)

spher.<th,ph> = A.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') 

spher_to_stereo = spher.transition_map(stereoN_A, (sin(th)*cos(ph)/(1-cos(th)) , sin(th)*sin(ph)/(1-cos(th))) )

spher_to_stereo.set_inverse(2*atan(1/sqrt(x^2+y^2)), atan2(y,x))
