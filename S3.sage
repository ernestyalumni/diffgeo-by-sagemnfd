## S3.sage
## This is my implementation of a 3-sphere, S^3, based upon 
## SM_sphere_S2 Sphere S^2 from Eric Gourgoulhon, Michal Bejger (2014) 
##   using these authors' sagemanifolds 0.7 package for Sage Math
## The main reference that I'll liberally copy from is the pdf SM_sphere_S2 on the authors' webpage 
##   for sagemanifolds
################################################################################################              
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                                                      
##                                                                                                            
## 20150413
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

## 20150413 EY : This is v2, version 2, because I am implementing a solution to 
## Tutorial 11 Symmetry, Exercise 3: Lie derivative, 
##   WE Heraeus International Winter School on Gravity and Light
##

S3 = Manifold(3, 'S^3', r'\mathbb{S}^3', start_index=1) ; print S3

#####
## Charts on S^3
#####

U = S3.open_subset('U') ; print U
V = S3.open_subset('V') ; print V

S3.declare_union(U,V) # declare U,V cover S^3

stereoN.<x,y,z>    = U.chart()
S = U.point((0,0,0), name='S') ; print S

stereoS.<xp,yp,zp> = V.chart("xp:x' yp:y' zp:z'")
N = V.point((0,0,0), name='N') ; print N

trans = stereoN.transition_map(stereoS, (x/(x^2+y^2+z^2), y/(x^2+y^2+z^2), z/(x^2+y^2+z^2)), intersection_name='W', restrictions1=x^2+y^2+z^2!=0, restrictions2=xp^2+yp^2+zp^2!=0)

W = U.intersection(V)

inv_trans = trans.inverse()

stereoN_W = stereoN.restrict(W)
stereoS_W = stereoS.restrict(W)

#####
## Scalar fields
#####

f0 = S3.scalar_field({stereoN: function('F',x,y)}, name='f_0') ; f0.display()
f0.add_expr_by_continuation(stereoS, W) ; f0.display()

#####
## Vector fields
#####

u = U.vector_field('u')
u[:] = [function('u_x',x,y,z), function('u_y', x,y,z) , function('u_z',x,y,z) ] ; u.display()

#####
## Embedding of S^3 into R^4
#####

R4 = Manifold(4, 'R^4', r'\mathbb{R}^4', start_index=0)
cart.<T,X,Y,Z> = R4.chart() ; cart

Phi = S3.diff_mapping(R4, {(stereoN, cart): [(x^2+y^2+z^2-1)/(1+x^2+y^2+z^2), 2*x/(1+x^2+y^2+z^2),2*y/(1+x^2+y^2+z^2),2*z/(1+x^2+y^2+z^2) ], 
    (stereoS, cart): [(1-xp^2-yp^2-zp^2)/(1+xp^2+yp^2+zp^2), 2*xp/(1+xp^2+yp^2+zp^2),2*yp/(1+xp^2+yp^2+zp^2),2*zp/(1+xp^2+yp^2+zp^2) ]}, name='Phi',latex_name=r'\Phi') ; Phi.display()

#####
## Riemannian metric on S^3
#####

h = R4.riemann_metric('h')
h[0,0], h[1,1], h[2,2], h[3,3] = 1,1,1,1 ; h.display()

g = S3.riemann_metric('g')
g.set( Phi.pullback(h) ) ; print g


#####
## Connection
#####
"""
nabla = g.connection()
# nabla.coef(spher.frame())[:,spher]
"""


#####
## Curvature
#####
"""
Riem = g.riemann() ; print Riem
Riem.display()
Riem.display( spher.frame(), spher)

Ric = g.ricci() ; Ric.display()

R = g.ricci_scalar() ; R.display()
"""

#####
## Spherical coordinates
#####

A = W.open_subset('A', coord_def={stereoN_W: (z!=0, y<0, x<0)})
stereoN_A = stereoN_W.restrict(A)

spher.<th,ph,ph3> = A.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi ph3:(0,2*pi):\phi^3') # cf. http://en.wikipedia.org/wiki/N-sphere

spher_to_stereo = spher.transition_map(stereoN_A, (sin(th)*sin(ph)*cos(ph3)/(1-cos(th)) , sin(th)*sin(ph)*sin(ph3)/(1-cos(th)) , sin(th)*cos(ph)/(1-cos(th) ) ) ) 

spher_to_stereo.set_inverse(2*atan(1/sqrt(x^2+y^2+z^2)), atan( sqrt(x^2+y^2)/z ) , atan2(y,x))

