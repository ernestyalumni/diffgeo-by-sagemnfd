## Rd.sage
## Euclidean Rd
## I've used and relied upon heavily the SM_tutorial.pdf SageManifolds tutorial by 
## Eric Gourgoulhon, Michal Bejger (2015)
## who both not only made a great Sage package and lucid documentation but are very helpful,
## and generous with helping anyone, and we're all richer for it; we really are.    
## 
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
##                                                            
## 20150826
##                                                                               
## This program, along with all its code, is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation; either version 2 of the License, or   
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
## Donate, and support my other scientific and engineering endeavors at 
## ernestyalumni.tilt.com                                                      
##                          
## Facebook     : ernestyalumni                                                   
## linkedin     : ernestyalumni                                                    
## Tilt/Open    : ernestyalumni                                                    
## twitter      : ernestyalumni                                                   
## youtube      : ernestyalumni                                                   
## wordpress    : ernestyalumni                                                    
## 
############################################################################

class Rd(object):
      """
      class Rd
      """
      def __init__(self,d):
      	  self.M = Manifold(d,'R'+str(d),r'\mathbb{R}^d',start_index=1)
	  self.X = self.M.chart(r" ".join([r"x"+str(i) for i in range(1,d+1)]))
	  self.U = self.M.open_subset('U',coord_def={self.X:[self.X[i] != 0 for i in range(1,d+1)]})
	  self.X_U = self.X.restrict(self.U)
	  if d==2:
	     self.sph = self.U.chart(r'rh:(0,+oo):\rho ph:(0,2*pi):\phi')
	  elif d==3:
	     self.sph = self.U.chart(r'rh:(0,+oo):\rho th:(0,pi):\theta ph:(0,2*pi):\phi')
	     self.cyl = self.U.chart(r'r:(0,+oo):r ph:(0,2*pi):\phi z')
	  elif d>3:
	     self.sph = self.U.chart(r'rh:(0,+oo):\rho '+r" ".join([r"th"+str(i)+r":(0,pi)" for i in range(1,d+1-2)])+r' ph:(0,2*pi):\phi')


	  if d==2:	  
	     self.transit_sph_to_X = self.sph.transition_map(self.X_U, [self.sph[1]*cos(self.sph[2]), self.sph[1]*sin(self.sph[2])])
	  elif d==3:
	     self.transit_sph_to_X = self.sph.transition_map(self.X_U, [self.sph[1]*sin(self.sph[2])*cos(self.sph[3]),self.sph[1]*sin(self.sph[2])*sin(self.sph[3]),self.sph[1]*cos(self.sph[2])])
	     self.transit_cyl_to_X = self.cyl.transition_map(self.X_U, [self.cyl[1]*cos(self.cyl[2]),self.cyl[1]*sin(self.cyl[2]),self.cyl[3]])
	  elif d>3:
	     self.transit_sph_to_X = self.sph.transition_map(self.X_U, [self.sph[1]*prod([sin(self.sph[i]) for i in range(2,j)])*cos(self.sph[j]) for j in range(3,d+1)]+[self.sph[1]*prod([sin(self.sph[i]) for i in range(2,d+1)]), self.sph[1]*cos(self.sph[2])])
	     
	  if d==2:
	     self.transit_sph_to_X.set_inverse(sqrt( self.X_U[1]^2+self.X_U[2]^2),atan2(self.X_U[2],self.X_U[1]))
	  elif d==3:
	     self.transit_sph_to_X.set_inverse(sqrt(self.X_U[1]^2+self.X_U[2]^2+self.X_U[3]^2),atan2(sqrt(self.X_U[1]^2+self.X_U[2]^2),self.X_U[3]),atan2(self.X_U[2],self.X_U[1]))
	     self.transit_cyl_to_X.set_inverse(sqrt(self.X_U[1]^2+self.X_U[2]^2),atan2(self.X_U[2],self.X_U[1]),self.X_U[3])
	  elif d>3:
	     gen_transit_list = [ sqrt(sum([ self.X_U[i]**2 for i in range(1,d+1)])), atan2( sqrt( sum([self.X_U[i]**2 for i in range(1,d)])), self.X_U[d]) ] +  [atan2( sqrt( sum( [ self.X_U[i]**2 for i in range(d-1,j,-1)])), self.X_U[j]) for j in range(1,d-2)] + [ 2*atan2( self.X_U[d-1] , self.X_U[d-2] + sqrt( (self.X_U[d-1])**2 + (self.X_U[d-2])**2)) ] 
	     self.transit_sph_to_X.set_inverse( *gen_transit_list)

      def equip_metric(self):
      	  self.g = self.M.riemann_metric('g')
	  for i in self.M.index_generator(1):
	      self.g[i[0],i[0]] = 1      	  

def make_pt(M):
    """
    make_pt = make_pt(M)
    INPUT
    M = sagemanifold Manifold

    EXAMPLES of USAGE
    p = make_pt(R3.M)
    """
    coords = M.default_chart()[:]
    farglst = ['p',]+list(coords)
    p = M.scalar_field( function(*farglst) )
    return p


def make_ut(M):
    """
    make_u = make_u(M)
    
    INPUT
    M = sage manifold Manifold
    """				  
    t = var('t')
    coords = M.default_chart()[:]
    # ucomplst components list of u  				  
    ucomplst = []
    for i in M.index_generator(1):
    	farglst = ['u'+str(i[0]),] + list(coords) + [t,]
	ucomplst.append( function( *farglst ) )
    u = M.vector_field(name='u')
    u[:] = ucomplst
    return u

def make_material_der(u, M):
    """
    make_material_der = make_material_der(u,M)

    EXAMPLES of USAGE:
    R3 = Rd(3)
    u3t = make_ut(R3.M)
    udu = make_material_der(u3t, R3.M)
    """
    uedcomp = []
    for ui in u.components()[:]:
	uedcomp.append( u( M.scalar_field( ui ) ))
    X = sum( [ uedcomp[i[0]-1]*M.default_frame()[i[0]] for i in M.index_generator(1) ] )
    return X


def div(u,g):
    """	
    div = div(u,g)
    Return the divergence of vector field u \in \mathfrak{X}(M), given the metric g for the manifold M
    """
    uflat = g['_ij']*u['^j']
    return xder( uflat.hodge_star(g) )

def grad(p,g):
    """
    grad = grad(p,g)

    EXAMPLE of USAGE
    R3 = Rd(3)
    p = make_pt(R3.M)
    grad(p,R3.g)
    """
    dp = xder( p )
    gradp = g.inverse()['^ij']*dp['_j']
    return gradp

def curl(u,g):
    """
    curl = curl(u,g)
    Return the curl of vector field u \in \mathfrak{X}(M), given the metric g for the manifold M
    """
    uflat = g['_ij']*u['^j']
    duflat = xder( uflat )
    return duflat.hodge_star(g)


##############################
## Usage Examples
##############################
"""
R3 = Rd(3)
R3.equip_metric()
p = make_pt(R3.M)
u3t = make_ut(R3.M)
R3.g(u3t,u3t)
R3.g(u3t,u3t).lie_der(u3t).display()

# 2 ways to contract tensors, or 2 ways to flat and sharp, through cotangent-tangent isomorphism:
R3.g.contract(0,u3t,0)
R3.g['_ij']*u3t['^j']

# then do d(u^{\flat})
duflat = xder(R3.g['_ij']*u3t['^j'])

# These steps are all wrapped up on curl = curl(u,g)
vorticity = curl( u3t, R3.g)
"""
