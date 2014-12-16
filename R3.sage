# R3.sage
# 3-dimensional Euclidean space, implementing Sage Manifold tutorial code examples 
# 3-dimensional Euclidean space, implemented using Sage Manifold tutorial code examples that is ready to load in Sage Math already. Switch to the directory and then run it.
#
# 20141216
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

R3 = Manifold(3,'R3',r'\mathbb{R}^3',start_index=1)
R3ch.<x,y,z> = R3.chart()

U = R3.open_domain('U',coord_def={R3ch: (x!=0,y!=0,z!=0)})
R3ch_U = R3ch.restrict(U)
sph.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi') 

transit_sph_to_R3ch = sph.transition_map(R3ch_U,[r*sin(th)*cos(ph),r*sin(th)*sin(ph),r*cos(th)])
transit_R3ch_to_sph = transit_sph_to_R3ch.set_inverse(sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y,x))


# vector frame on this chart
R3ch.frame()
# X.default_frame()
E = R3ch.frame()

#
# The Metric : Minkowski flat metric
#

g = R3.riemann_metric('g')
print g

# defining the metric
g[1,1], g[2,2], g[3,3] = 1, 1, 1
g.view()

nabla = g.connection()
print nabla ; nabla
nabla.coef()[:]


#
# volume form
#

vol_3 = g.volume_form()
vol_3.view()


#
# differential forms
# 



#
# cf. Sh. Sternberg Semi-Riemannian Geometry and General Relativity
# Ch. 11 Star
# 11.3.2 for R^3
# pp. 241
#

g.volume_form().view() # eps_g = dx/\dy/\dz # EY !!! reproduces the volume form !!!
g.volume_form().hodge_star(g).view()
g.volume_form().hodge_star(g).hodge_star(g).view()

R3ch.coframe()[1].hodge_star(g).view()
R3ch.coframe()[2].hodge_star(g).view()
R3ch.coframe()[3].hodge_star(g).view()
Theta = R3ch.coframe()

a = function('a',x,y,z)
b = function('b',x,y,z)
c = function('c',x,y,z)
A = function('A',x,y,z)
B = function('B',x,y,z)
C = function('C',x,y,z)

theta_110302 = a * Theta[1] + b * Theta[2] + c * Theta[3] 
omega_110302 = A * Theta[3].hodge_star(g)  - B * Theta[2].hodge_star(g) + C * Theta[1].hodge_star(g)

# taking a look at each piece of what we want to show for Exercise (Problem) 8 on pp. 241, (11.10), d*\theta is interesting because it's the divergence, essentially, and taking the star again makes it the divergence, for real

theta_110302.hodge_star(g).view()
xder( theta_110302.hodge_star(g) ).view()

# d * d * \theta piece is interesting as it contains some "mixed" partial derivatives
xder( xder( theta_110302.hodge_star(g) ).hodge_star(g) )

# for (11.10) *d*d\theta - d*d*\theta

Eq11010 = xder( xder( theta_110302 ).hodge_star(g) ).hodge_star(g) - xder( xder( theta_110302.hodge_star(g) ).hodge_star(g)) # Equation 11.10

print Eq11010.view() # this is what we want, it's like the Laplacian is applied to each component

Eq11011 = - xder( xder( omega_110302 ).hodge_star(g) ).hodge_star(g) + xder( xder( omega_110302.hodge_star(g) ).hodge_star(g) )

print Eq11011.view()

#
# the following code is simply to demonstrate how powerful, I believe, sagemanifolds is, in its symbolic computation power: what I see on "paper," can be computed out and a lot of tedious algebra is avoided, and it makes differential geometry "fun", and fun is most important to me (if math and physics stopped being fun, then I wouldn't do it)
# 

print "\n Let's reproduce the volume form"
print g.volume_form().view()

print "\n Let's take the Hodge star of each of the basis covectors or 1-forms"
print R3ch.coframe()[1].hodge_star(g).view()
print R3ch.coframe()[2].hodge_star(g).view()
print R3ch.coframe()[3].hodge_star(g).view()

print "\n We're given theta 1-form and Omega 2-form; let's reproduce that:"
print theta_110302.view()
print omega_110302.view()

print "\n For Exercise 8 (Problem 8) in 11.3.2 of S. Sternberg's Semi-Riemannian Geometry and General Relativity is easily done, without resorting to tedious calculus"
print "\n This is the LHS of Equation (11.10)"
print "xder( xder( theta_110302 ).hodge_star(g) ).hodge_star(g) - xder( xder( theta_110302.hodge_star(g) ).hodge_star(g))"
print "\n This is the RHS of Equation (11.10)"
print Eq11010.view()

print "\n This is the LHS of Equation (11.11)"
print "- xder( xder( omega_110302 ).hodge_star(g) ).hodge_star(g) + xder( xder( omega_110302.hodge_star(g) ).hodge_star(g) )"
print "\n This is the RHS of Equation (11.11)"
print Eq11011.view()
