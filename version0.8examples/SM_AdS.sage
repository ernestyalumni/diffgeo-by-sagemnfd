## SM_AdS.sage
## This is my implementation of Anti-de Sitter spacetime
##   utilizing
## Sage Math
##  and sagemanifolds package 
## Much if not all credit should go to Eric Gourgoulhon, Michal Bejger for solely creating an 
## incredible and mathematically sound Sage package; and on top of that, they genuinely help others,## without knowing who you are, and those guys are really worth their weight in gold; because 
## the fact is that they aren't that many people like that in the world, and we need them. 
##
## The main reference that I'll liberally copy from is from is 
##   SM_AdS.pdf
##   http://sagemanifolds.obspm.fr/examples/pdf/SM_AdS.pdf
################################################################################
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>
## 
## 20150808
## 
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, 
## or (at your option) any later version. 
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
## Governing the ethics of using this program, I default to the Caltech 
## Honor Code:                           
## ``No member of the Caltech community shall take unfair advantage of 
## any other member of the Caltech community.''
## 
## Donate and support my scientific and engineering efforts here:
## ernestyalumni.tilt.com 
## 
## Facebook     : ernestyalumni 
## linkedin     : ernestyalumni 
## Tilt/Open    : ernestyalumni 
## twitter      : ernestyalumni 
## youtube      : ernestyalumni 
## wordpress    : ernestyalumni 
##  
################################################################################

#########################
## Spacetime manifold
#########################

M = Manifold(4, 'M', r'\mathcal{M}')
print M; M

M0 = M.open_subset('M_0', r'\mathcal{M}_0')
X_hyp.<ta,rh,th,ph> = M0.chart(r'ta:\tau rh:(0,+oo):\rho th:(0,pi):\theta ph:(0,2*pi):\phi')
print X_hyp; X_hyp

###################################
## \mathbb{R}^5 as an ambient space
###################################

R5 = Manifold(5,'R5',r'\mathbb{R}^5')
X5.<U,V,X,Y,Z> = R5.chart()
print X5; X5

var('b')
assume(b>0)
Phi = M.diff_mapping(R5, [sin(b*ta)/b*cosh(rh), 
      			 cos(b*ta)/b * cosh(rh),
			 sinh(rh)/b*sin(th)*cos(ph),
			 sinh(rh)/b*sin(th)*sin(ph),
			 sinh(rh)/b*cos(th)],
			 name='Phi', latex_name=r'\Phi')
print Phi; Phi.display()

p = M.point((ta,rh,th,ph),name='p') ; print p
p.coord()

q = Phi(p); print q
q.coord()

(Uq,Vq,Xq,Yq,Zq) = q.coord()
s = -Uq^2 - Vq^2 + Xq^2 + Yq^2 + Zq^2
s.simplify_full()

graph1 = X_hyp.plot(X5, mapping=Phi, ambient_coords=(V,X,U), fixed_coords={th:pi/2,
ph:0}, ranges={ta:(0,2*pi), rh:(0,2)}, nb_values=9, color={ta:'red', rh:'brown'},
thickness=2, parameters={b:1}, label_axes=False)
graph2 = X_hyp.plot(X5, mapping=Phi, ambient_coords=(V,X,U), fixed_coords={th:pi/2,
ph:pi}, ranges={ta:(0,2*pi), rh:(0,2)}, nb_values=9, color={ta:'green', rh:'brown'},
thickness=2, parameters={b:1}, label_axes=False)
show(set_axes_labels(graph1+graph2,'V','X','U'), aspect_ratio=1)

####################
## Spacetime metric
####################

h = R5.metric('h')
h[0,0], h[1,1], h[2,2], h[3,3], h[4,4] = -1, -1, 1, 1, 1
h.display()

g = M.metric('g')
g.set( Phi.pullback(h) )

g.display()

g[:]

############### 
## Curvature
############### 

Riem = g.riemann()
print Riem
Riem.display()

Riem.display_comp(only_nonredundant=True)

Ric = g.ricci()
print Ric
Ric.display()

Ric[:]

R = g.ricci_scalar()
print R
R.display()

delta = M.tangent_identity_field()
Riem == - (R/6)*(g*delta).antisymmetrize(2,3) # 2,3 = last positions of the type-(1,3) tensor g*delta

Lambda = -3*b^2
Ric - 1/2*R*g + Lambda*g == 0 

#########################
## Spherical coordinates
#########################

X_spher.<ta,r,th,ph> = M0.chart(r'ta:\tau r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
print X_spher ; X_spher

hyp_to_spher = X_hyp.transition_map(X_spher, [ta,sinh(rh)/b, th,ph])
hyp_to_spher.display()

hyp_to_spher.set_inverse(ta, asinh(b*r), th,ph)
spher_to_hyp = hyp_to_spher.inverse()
spher_to_hyp.display()

g.display(X_spher.frame(), X_spher)

