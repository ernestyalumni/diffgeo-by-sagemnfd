# Min4d.sage
# 4-dimensional Minkowski space, implementing Sage Manifold tutorial code examples 
# 4-dimensional Minkowski space, implemented using Sage Manifold tutorial code examples that is ready to load in Sage Math already. Switch to the directory and then run it.
#
# This sagemanifolds code helps solve these book problems:
# Wald 2.8b+
# (note that previous problem is from Problem Set 3 of Ph236a, taught by Kamionkowski at Caltech, Oct. 17, 2006
#
# 20141005
# 
#
# Fund Science! & Help Ernest finish his Physics Research! : quantum super-A-polynomials - Ernest Yeung
#                                               
# http://igg.me/at/ernestyalumni2014                                                                             
#                                                              
# Facebook     : ernestyalumni  
# github       : ernestyalumni                                                                     
# gmail        : ernestyalumni                                                                     
# linkedin     : ernestyalumni                                                                             
# tumblr       : ernestyalumni                                                               
# twitter      : ernestyalumni                                                             
# youtube      : ernestyalumni                                                                
# indiegogo    : ernestyalumni                                                                        
# 
# Ernest Yeung was supported by Mr. and Mrs. C.W. Yeung, Prof. Robert A. Rosenstone, Michael Drown, Arvid Kingl, Mr. and Mrs. Valerie Cheng, and the Foundation for Polish Sciences, Warsaw University.  
#
# This code is open-source, governed by the Creative Common license.  Use of code is governed by the Caltech Honor Code: ``No member of the Caltech community shall take unfair advantage of any other member of the Caltech community.'' 
#
# NOTE : this sage code is for Sage 6.3. and sagemanifolds 0.5
# if you need help upgrading, don't hesitate to ask
# 
# NOTE : I'm going to take the liberty to quote word for word, directly from the sage manifolds tutorial
# 

R4 = Manifold(4,'R4',r'\mathbb{R}^4',start_index=0)
flch.<t,x,y,z> = R4.chart() # "flat chart"

U = R4.open_domain('U',coord_def={X: (x!=0,y!=0,z!=0)})
flch_U = flch.restrict(U)
sph4.<t,r,th,ph> = U.chart(r't:(-oo,oo) r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi') 
cyl4.<t,r,ph,z>  = U.chart(r't:(-oo,oo) r:(0,+oo) ph:(0,2*pi):\phi z:(-oo,oo)')

transit_sph4_to_flch = sph4.transition_map(flch_U,[t,r*sin(th)*cos(ph),r*sin(th)*sin(ph),r*cos(th)])
transit_flch_to_sph4 = transit_sph4_to_flch.set_inverse(t,sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y,x))

transit_cyl4_to_flch = cyl4.transition_map(flch_U,[t,r*cos(ph), r*sin(ph), z] )
transit_flch_to_cyl4 = transit_cyl4_to_flch.set_inverse(t,sqrt(x^2+y^2), atan2( y,x), z)


# vector frame on this chart
flch.frame()
# X.default_frame()
E = flch.frame()

#
# The Metric : Minkowski flat metric
#

g = R4.riemann_metric('g')
print g

# defining the metric
g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
g.view()
g.view(sph4.frame(), sph4) # g = -dt*dt + dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
g.view(cyl4.frame(), cyl4) # g = -dt*dt + dr*dr + r^2 dph*dph + dz*dz # EY : 20141221 !!! This is great stuff


nabla = g.connection()
print nabla ; nabla
nabla.coef()[:]


#
# volume form
#

vol_4 = g.volume_form()
vol_4.view()


#
# differential forms
# 


#
# The following problem is from Ph236a taught by Kamionkowski at Caltech, Problem Set 3, Oct. 17, 2006 Problems 6
# 6. (Wald 2.8b+) Rotating coordinates
# 

# cf. http://www.sagemath.org/doc/reference/constants/sage/symbolic/constants.html
from sage.symbolic.constants import Constant
omega = var('omega')

# cf. http://www.sagemath.org/doc/reference/calculus/sage/calculus/calculus.html
F(t,x,y,z) = [ t , (x**2 + y**2)**(1/2) * cos( arctan2( y,x) - omega *t ) , (x**2 + y**2)**(1/2) * sin( arctan2( y,x) - omega * t) , z ]

F.diff()
diff(F) # print this out and it looks like a mess; use full_simplify on each row



print "\n Ph236a 2006-2007 Caltech Kamiokowski" 
print "P.S.3, Problem 6, (Wald 2.8b+) Rotating coordinates"
print "6.(a)"

print "\n This is DF"
print diff(F)(t,x,y,z)

print "\n It looks like a mess; in Sage Math/Python, use full_simplify"
for i in range(0,4):
    print [ j.full_simplify() for j in diff(F).rows()[i](t,x,y,z) ]

print "\n Take the inverse of this DF and do full_simplify"
for i in range(0,4):
    print [ j.full_simplify() for j in diff(F).inverse().rows()[i](t,x,y,z) ]


gprime = (  (matrix(g[:]) * diff(F).inverse()).T * diff(F).inverse() )(t,x,y,z)

print "\n gprime is the following:"
for i in range(0,4):
    print [j.full_simplify() for j in gprime.rows()[i] ]


print "\n 6.(b)"
Fb(t,r,ph,z) = [ t , r , ph - omega *t , z]
gbprime = ( ( diagonal_matrix( [ -1 , 1 , r^2, 1] ) * diff(Fb).inverse() ).T * diff(Fb).inverse() )(t,r,ph, z)

print "\n the Riemannian metric, or 'line element' (EY : 20141221 in my opinion, 'line element' is a misnomer) is the following:"
print gbprime
