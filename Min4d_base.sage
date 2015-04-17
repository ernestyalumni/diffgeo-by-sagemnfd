## Min4d_base.sage
## 4-dimensional Minkowski space, implementing Sage Manifold 0.7 tutorial code examples 
## 4-dimensional Minkowski space, implemented using Sage Manifold 0.7 tutorial code examples 
## that is ready to load in Sage Math already. Switch to the directory and then run it.
## base is for just the base: most essentially functions
## 
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

R4 = Manifold(4,'R4',r'\mathbb{R}^4',start_index=0)
cart.<t,x,y,z> = R4.chart() # "flat chart"


U = R4.open_subset('U',coord_def={cart: [x!=0,y!=0,z!=0]})  # [x!=0 AND y!=0 AND z!=0]

cart_U = cart.restrict(U)
sph.<t,r,ph,th> = U.chart(r't:(-oo,oo) r:(0,+oo) ph:(0,2*pi):\phi th:(0,pi):\th') 
cyl.<t,r,ph,z>  = U.chart(r't:(-oo,oo) r:(0,+oo) ph:(0,2*pi):\phi z:(-oo,oo)')

transit_sph_to_cart = sph.transition_map(cart_U,[t,r*sin(th)*cos(ph),r*sin(th)*sin(ph),r*cos(th)])
transit_cart_to_sph = transit_sph_to_cart.set_inverse(t,sqrt(x^2+y^2+z^2), atan2(y,x), atan2(sqrt(x^2+y^2),z))

transit_cyl_to_cart = cyl.transition_map(cart_U,[t,r*cos(ph), r*sin(ph), z] )
transit_cart_to_cyl = transit_cyl_to_cart.set_inverse(t,sqrt(x^2+y^2), atan2( y,x), z)


# vector frame on this chart
ecart = cart.frame()

#
# The Metric : Minkowski flat metric
#

g = R4.lorentz_metric('g')
print g

# defining the metric
g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
g.display()
g.display(sph.frame(), sph) # g = -dt*dt + dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
g.display(cyl.frame(), cyl) # g = -dt*dt + dr*dr + r^2 dph*dph + dz*dz # EY : 20141221 !!! This is great stuff

g_U = U.lorentz_metric('g_U')
g_U[0,0], g_U[1,1], g_U[2,2], g_U[3,3] = -1, 1, 1, 1

nabla = g.connection()
# print nabla ; nabla
# nabla.coef()[:]


#
# volume form
#

vol_4 = g.volume_form()
vol_4.display()


#
# differential forms
# 


