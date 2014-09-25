# EY_sm_tut_wkbk.sage
# Sage Manifold tutorial code examples 
# Sage Manifold tutorial code examples that is ready to load in Sage Math already. Switch to the directory and then run it.
#
# 20140924
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
# Defining a manifold
#
M = Manifold(3, 'M', r'\mathcal{M}', start_index=1)
M
print M
type(M)
for i in M.irange():
    print i

M0 = Manifold(3,'M', r'\mathcal{M}')
for i in M0.irange():
    print i

#
# Defining a chart on the manifold
#

X.<x,y,z> = M.chart()
X[1]
X[2]
X[3]
X[:]
z is X[3]

# Functions of the chart coordinates

f = X.function(x+y^2+z^3) ; f
f.view()

f(1,2,3)

type(f)

f0(x,y,z) = cos(x)^2 ; g0(x,y,z) = sin(x)^2

f0 + g0

f1 = X.function(cos(x)^2) ; g1 = X.function(sin(x)^2)
f1 + g1

# EY 20140924
f1 is f0 # False, EY : they're different entities
f1 == f0 # True, same output, they do the same thing

(f0 + g0).simplify_trig()

f0
f1
f1.view()
f1.expr()
type(f1.expr())

#
# Introducing a second chart on the manifold
#

U = M.open_domain('U', coord_def={X: (y!=0, x<0)})

X_U = X.restrict(U); X_U

Y.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi'); Y

th, ph

Y[2], Y[3]

assumptions()
simplify(abs(r))

simplify(abs(r))
simplify(abs(x)) # no simplification occurs since x can take any value in R

transit_Y_to_X = Y.transition_map(X_U, [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th) ])
transit_Y_to_X

transit_Y_to_X(r,th,ph)

transit_X_to_Y = transit_Y_to_X.set_inverse(sqrt(x^2+y^2+z^2),atan2(sqrt(x^2+y^2),z),atan2(y,x))

M.atlas()
M.default_chart()

U.atlas()
U.default_chart()

#
# Points on the manifold
#
p = M.point((1,2,-1),X,name='p') ; print p ; p

# Since X = (M,(x,y,z)) is the manifold's default chart, its name can be omitted
p = M.point((1,2,-1), name='p') ; print p ; p

p in M

p in U

p.coord()

q = M.point((1,0,2), name='q')

q in U
q.coord(Y)

p.coord(Y)
q==p
p1=U.point((sqrt(6),pi-atan(sqrt(5)),atan(2)),Y)

p.parent()
q.parent()
p1.parent()

#
# Scalar fields
#

EY : keep in mind that so far, U has restrictions, 

f = U.scalar_field({X_U: x+y^2+z^3}, name='f') ; print f
# since in the present case, there's only 1 chart in the dictionary, an alternative writing is the following:
f = U.scalar_field(x+y^2+z^3, chart=X_U, name='f'); print f 
# Since X_U is the domain's default chart, it can be omitted in the above declaration
f = U.scalar_field(x+y^2+z^3, name='f') ; print f

f(p)
f.view(X_U)
f.view()
f.view(Y)

f.function_chart(X_U)
f.function_chart(X_U).view()
f.function_chart(Y)
f.function_chart(Y).view()

f.expr(X_U)
f.expr(Y)
f.expr(Y) is f.function_chart(Y).expr()

# A scalar field can also be defined by some unspecified function of the coordinates:
h = U.scalar_field(function('H', x, y, z), name='h') ; print h

h.view()
h.view(Y)
h(p)

# The parent of f is the set C^{\infty}(U) of all smooth scalar fields on U, which is a commutative algebra over \mathbb{R}
# This is exciting because now we can do the C^{\infty}(U) algebra in a civilized and sane manner
CU = f.parent() ; CU
print CU
CU.category()
CU.base_ring()

s = f + 2*h ; print s

s.view()

#
# Tangent spaces
# 

Tp = p.tangent_space() ; Tp
print Tp
Tp.category()
Tp.dim()

Tp.bases()

Tq = q.tangent_space()
Tq.bases()

# EY : a random element? I think the appropriate wording is an arbitrary element, this is simply semantics
v = Tp.an_element() ; print v

v.view()

u = Tq.an_element() ; print u
u.view()
v.parent()

u.parent()
s = u + v

#
# Vector Fields
#

# EY : remember, X was a chart on M
X.frame()
X.frame().domain()

Y.frame()
Y.frame().domain()

M.frames()
U.frames()
M.default_frame()
M.default_frame() is M.default_chart().frame()
U.default_frame() is U.default_chart().frame()

e = U.default_frame() ; e2 = e[2] ; e2
print e2
v = e[2] + 2*x*e[3] ; print v
v.view()

v = U.vector_field(name='v') # vector field defined on the open domain U
v[1] = 1+y
v[2] = -x
v[3] = x*y*z
v.view()

v.parent()
print v.parent()

print v.parent().category()

v.parent().base_ring()

f.view()

# vector field acts on scalar fields
s = v(f) ; print s  

s.view()

e[3].view()
e[3](f).view()  # EY : this is hot, a vector does on f, C^{\infty} function what it should
w = U.vector_field(name='w')
w[2] = 3 
w.view() # unset components are assumed to be zero
v.view(Y.frame())

# one should add the corresponding chart as the second argument
v.view(Y.frame(),Y)
for i in M.irange(): e[i].view(Y.frame(),Y)

v[:]
v.comp(Y.frame())[:]
v[Y.frame(),:] # shortcut
v[Y.frame(),:,Y]

print v[[1]]
v[[1]].view()

v[[1]].expr(X_U)

# A vector field can be defined with components being unspecified functions of the coordinates:
u = U.vector_field(name='u')
u[:] = [function('u_x', x,y,z), function('u_y', x,y,z), function('u_z', x,y,z)]
u.view()

# Values of vector fields at a given point

vp = v.at(p) ; print vp

vp.view()
p.coord(X_U)
v.view(X_U.frame(),X_U)
vp.set_name(latex_name='v|_p')
vp.view()

vp.parent()
vp in p.tangent_space()

up = u.at(p) ; print up
up.view()

#
# 1-forms
#

df = f.differential() ; print df
df.view()

dX = X.coframe() ; dX

M.coframes()

dfp = df.at(p) ; print dfp
dfp.view()

p.coord()
dfp.parent()
dfp.parent() is p.tangent_space().dual()

print vp ; vp.view()

dfp(vp)

print up ; up.view()
dfp(up)

h.view() ; dh = h.differential() ; dh.view()

# modified Pynac output
from sage.symbolic.pynac import textbook_style_deriv, omit_function_args
textbook_style_deriv(True)
omit_function_args(True)

dh.view()

# A 1-form can also be defined from scratch:
om = U.one_form('omega', r'\omega') ; print om
om[:] = [x^2+y^2, z, x-z] # components in the default coframe (dx,dy,dz)
om.view()

om[Y.frame(), :, Y] = [r*sin(th)*cos(ph), 0, r*sin(th)*sin(ph)]
om.view(Y.frame(), Y)
om.view()
om[:] = [x^2+y^2, z, x-z]
om.view()
om.view(Y.frame(),Y)

v.view() ; om.view() ; print om(v) ; om(v).view()

df.view() ; print df(v) ; df(v).view()

u.view() ; om(u).view()

df.parent()
print df.parent()
print om.parent()

#
# Differential forms and exterior calculus
#

a = om.wedge(df) ; print a ; a.view()
a[:]
a.view(Y.frame(),Y)

a.set_name('A')
print a(i,v) ; a(u,v).view()

a(u,v) == - a(v,u)

a.symmetries()

# The exterior derivative of a differential form is taken by invoking the function xder():
dom = xder(om) ; print dom ; dom.view()

da = xder(a) ; print da ; da.view()

ddf = xder(df) ; ddf.view()
ddom = xder(dom) ; ddom.view()

#
# Lie derivative
# 

lv_om = om.lie_der(v) ; print lv_om ; lv_om.view()

lu_dh = dh.lie_der(u) ; print lu_dh ; lu_dh.view()

om.lie_der(v) == v.contract(xder(om)) + xder(om(v))

a.lie_der(v) == v.contract(xder(a)) + xder(v.contract(a)) # EY : amazing

v.lie_der(u)(f) == u(v(f)) - v(u(f)) # False

#
# Tensor fields of arbitrary rank
#

t = U.tensor_field(1,2, name='T') ; print t

t[1,2,1] = 1 + x^2
t[3,2,1] = x*y*z

t.view()

t[:]

print t[[1,2,1]] ; t[[1,2,1]].view()

print t[1,2,1] ; t[1,2,1]

print t(om, u, v) ; t(om, u, v).view()

t[Y.frame(), 1,1,1, Y]  # tensor components

#
# Tensor calculus
#

v.tensor_type() ; a.tensor_type()

# The tensor product is denoted by *
b = v*a ; print b ; b

a.symmetries()
b.symmetries()

s = - t + 2*f* b ; print s

c = t.self_contract(0,1)

tv = t.contract(2, v, 0)
print tv

tv1 = t.contract(v)
print tv1

tv1 == tv

#
# Metric structures
#

g = M.riemann_metric('g')
print g

g.parent()
print g.parent()
g.symmetries()
g[1,1], g[2,2], g[3,3] = 1,1,1
g.view()

g.view(Y.frame(), Y)

u.view() ; v.view(); print g(u,v) ; g(u,v).view()

# The Levi-Civita connection associated to the metric g:
nabla = g.connection()
print nabla ; nabla

# The Christoffel symbols with respect to the manifold's default coordinates:

nabla.coef()[:]

nabla.coef(Y.frame())[:, Y]

# The connection acting as a covariant derivative EY : !!! this is great stuff
nab_v = nabla(v)
print nab_v ; nab_v.view()

print nabla.torsion() ; nabla.torsion().view()

print nabla.riemann() ; nabla.riemann().view()

# consider a non-flat metric, change g_{rr}

g[Y.frame(), 1,1, Y] = 1/(1+r^2)
g.view(Y.frame(), Y)

# For convenience, change default chart on domain U to Y
U.set_default_chart(Y)

g.view(Y.frame())

g.view(X_U.frame(), X_U)

g[X_U.frame(), : , X_U]

g.add_comp_by_continuation(X.frame(),U,X)
g.view()

nabla = g.connection()

nabla.coef()[:]

nabla.coef(Y.frame())[:]

Riem = nabla.riemann()  # this took a while for my computer
print Riem ; Riem.view(Y.frame())

Ric = g.ricci()
print Ric ; Ric.view(Y.frame())

# Weyl tensor 
C = g.weyl() # this took a while for my computer to run, but not as long as Riem
print C ; C.view()

R = g.ricci_scalar()
print R ; R.view()

# Tensor transformations induced by a metric

print t
t.view()
t.view(X_U.frame(), X_U)

# Raising the last index of T with g
s = t.up(g, 2)
print s

s = t.up(g)
print s

s = t.down(g)
print s

#
# Hodge duality
# 

epsilon = g.volume_form()
print epsilon ; epsilon.view()

epsilon.view(Y.frame())

print f ; f.view()

sf = f.hodge_star(g)
print sf ; sf.view()

sf == f * epsilon.restrict(U)

print om ; om.view()

som = om.hodge_star(g)
print som ; som.view()

print a

sa = a.hodge_star(g)
print sa ; sa.view()

print da ; da.view()
sda = da.hodge_star(g)
print sda ; sda.view()

sf.hodge_star(g) == f
som.hodge_star(g) == om
sa.hodge_star(g) == a
sda.hodge_star(g) == da




