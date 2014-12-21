# R3.sage
# 3-dimensional Euclidean space, implementing Sage Manifold tutorial code examples 
# 3-dimensional Euclidean space, implemented using Sage Manifold tutorial code examples that is ready to load in Sage Math already. Switch to the directory and then run it.
#
# This sagemanifolds code helps solve these book problems:
# MTW, Problem 9.6, MTW is Misnor, Thorne, Wheeler, Gravitation
# MTW, Problems 9.7-8
# Wald 2.8b+
# (note that previous 3 problems are from Problem Set 3 of Ph236a, taught by Kamionkowski at Caltech, Oct. 17, 2006
# cf. Sh. Sternberg Semi-Riemannian Geometry and General Relativity,  Ch. 11 Star, 11.3.2 for R^3, pp. 241 Exercise 8 (Problem 8)
# Also a Levi-Civita symbol is implemented as a totally antisymmetric tensor 
# 
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
# The following 3 problems are from Ph236a taught by Kamionkowski at Caltech, Problem Set 3, Oct. 17, 2006 Problems 4-6
#
# 4. (MTW, Problem 9.6) Practice with dual bases
print "\n Ph236a 2006-2007 Caltech Kamiokowski" 
print "P.S.3, Problem 4, MTW (Misner Thorne Wheeler Gravitation) Problem 9.6, Practice with dual bases"
print "sph, which is what I called the chart on $\mathbb{R}^3$ gives the coordinate basis. sage/python code is sph.frame()[1].view() and plug in 2,3 for the $\theta$, and $\phi$ corresponding vector field \n "
for i in range(1,3):
      print sph.frame()[i].view()

# sagemanifolds note: EY : 20141217 I foudn that one should start from the open domain of a chart and then define the vector field, cf. Vector Fields section of SM_tutorial
ehatr  = U.frames()[1][1]
ehatth = (1/r) * U.frames()[1][2]
ehatph = (1/(r*sin(th) ) ) * U.frames()[1][3]
eorthon = [ ehatr , ehatth , ehatph ]

print "\n orthonormal basis: "
print ehatr.view(sph.frame())
print ehatth.view(sph.frame())
print ehatph.view(sph.frame())

orthon_1form_basis = [ U.coframes()[1][1] , 
		   r * U.coframes()[1][2] ,
		   r * sin(th) * U.coframes()[1][3] ]

print "\n This is the corresponding one-form basis $\{ \omega^{\widehat{r}}, \omega^{\widehat{\theta}}, \omega^{\widehat{ \phi } } \}$ "
for omega in orthon_1form_basis:
    print omega.view(sph.frame())

print "\n Evaluating each 1-form basis dual to the orthonormal on each orthonormal basis element (9 possibilities), gives us the following, printing out the identity matrix, essentially"
for omega in orthon_1form_basis:
    for e in eorthon:
    	print omega(e).view()

print "\n P.S.3, Problem 5 (MTW, Problems 9.7-8) Commutators for Euclidean space in spherical coordinates"
print "part (a)"
# define any ol' scalar field
fcar = U.scalar_field(function('fcar',x,y,z), name='fcar')
fsph = U.scalar_field(function('fsph', r,th,ph), name='fsph')

print "\n Here is the evaluation of the commutator, 1 being theta hat, 2 being phi hat"
print (  eorthon[1]( eorthon[2]( fsph ) ) - eorthon[2]( eorthon[1]( fsph ) ) ).view()

#
# making a Levi-Civita symbol
# EY : 20141218 Let's implement a Levi-Civita symbol as a totally antisymmetric rank-3 tensor
eps3 = R3.tensor_field(0,3,'epsilon', antisym=((0,1),(1,2),(0,2)) )
eps3[1,2,3] = 1

# indeed this gives the desired symbol:
[ [ [ ( ( i, j , k ) , eps3[i, j , k ] ) for i in range(1,4) ] for j in range(1,4) ] for k in range(1,4) ]

# Angular momentum operators
def L_i(i, cht=R3ch ): 
    """
    L_i angular momentum operator in ith component
    INPUTS : i is an integer, 1,2,3
    cht - chart, the default chart is the Cartesian chart for R3 Euclidean 3 dimensional space, R3ch
    
    """
    L_i_temp = sum( [ sum( [ eps3[ i , j , k ].expr() * R3ch[ j ] * R3ch.frame()[ k ] for k in range(1,4) ] ) for j in range(1,4) ] )
    return L_i_temp

print "\n The Angular momentum operators are the following (and note the implementation of a Levi-Civita symbol)"
for i in range(1,4):
    print L_i(i).view()

print "\n Commutators for the 3 'angular momentum operators' "
print ( L_i(1)( L_i(2)( fcar ) ) - L_i(2)( L_i(1)( fcar ) ) ).view()
print ( L_i(3)( L_i(1)( fcar ) ) - L_i(1)( L_i(3)( fcar ) ) ).view()
print ( L_i(2)( L_i(3)( fcar ) ) - L_i(3)( L_i(2)( fcar ) ) ).view()


#
# cf. Sh. Sternberg Semi-Riemannian Geometry and General Relativity
# Ch. 11 Star
# 11.3.2 for R^3
# pp. 241
#

print "\n cf. Sh. Sternberg Semi-Riemannian Geometry and General Relativity"

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
