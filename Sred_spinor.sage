## Sred_spinor.sage
## This is my implementation of Pauli matrices, spinors, spin helicity 
##   using Sage Math
## The main reference 
##   Mark Srednicki, Quantum Field Theory, Cambridge University Press, Sec. 50
## You'll need spinhelicity.sage to load in your current directory on Sage Math
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
## Sred_spinor.sage
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

load("spinhelicity.sage")

#####################################################################################
## EY : What follows is a demonstration of using spinhelicity.sage to solve the 
## problems in Sec. 50 Massless particles and spinor helicity in Mark Srednicki's
## Quantum Field Theory
#####################################################################################

sq_p_ket = make_square_ket( pptsph.coord(sph_k) )
bra_p_sq = make_square_bra( pptsph.coord(sph_k) )
an_p_ket = make_angle_ket( pptsph.coord(sph_k) )
bra_p_an = make_angle_bra( pptsph.coord(sph_k) )

sq_p_ket_cart = make_square_ket( ppt.coord(sph_k) )
bra_p_sq_cart = make_square_bra( ppt.coord(sph_k) )
an_p_ket_cart = make_angle_ket( ppt.coord(sph_k) )
bra_p_an_cart = make_angle_bra( ppt.coord(sph_k) )

sq_k_ket = make_square_ket( kptsph.coord(sph_k) )
bra_k_sq = make_square_bra( kptsph.coord(sph_k) )
an_k_ket = make_angle_ket( kptsph.coord(sph_k) )
bra_k_an = make_angle_bra( kptsph.coord(sph_k) )

sq_k_ket_cart = make_square_ket( kpt.coord(sph_k) )
bra_k_sq_cart = make_square_bra( kpt.coord(sph_k) )
an_k_ket_cart = make_angle_ket( kpt.coord(sph_k) )
bra_k_an_cart = make_angle_bra( kpt.coord(sph_k) )


#####
## Srednicki Prob 50.1
#####

bra_k_an * sq_p_ket # 0
bra_k_sq * an_p_ket # 0

RHS_50_35 = an_p_ket_cart * bra_p_sq_cart + sq_p_ket_cart * bra_p_an_cart
# RHS_50_35 = Matrix( [ [ RHS_50_35[i,j].simplify_full().trig_expand(half_angles="True") for i in range(0,4) ] for j in range(0,4) ] )
LHS_50_35 = -pslash( ppt.coord() )
LHS_50_35_v2 = -pslash( pptsph.coord( cart_k) ); LHS_50_35_v2 = LHS_50_35_v2.subs(pr==p0)

LHS_50_35[0:2,2:4] # row first, column next

# EY : 20150417
# commands I've tried:
# sageobj( RHS_50_35[0:2,2:4][0,1]._maxima_().demoivre() )
RHS_50_35[0:2,2:4][0,0].trig_expand(half_angles="True")
