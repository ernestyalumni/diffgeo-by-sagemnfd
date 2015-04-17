## SO3.sage
## This is my implementation of SO(3) Lie Group and so(3) Lie Algebra
## SM_tutorial from Eric Gourgoulhon, Michal Bejger (2015) 
##   using these authors' sagemanifolds 0.7 package for Sage Math
## The main reference that I'll liberally copy from is from SM_tutorial
##   on the authors' webpage 
##   for sagemanifolds
################################################################################
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>
## 
## 20150416
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
## If you like what I'm doing and would like to help and contribute support,
## please take a look at my crowdfunding campaign at ernestyalumni.tilt.com 
## and Patreon, read my mission statement and give your financial support, 
## no matter how small or large, if you can and to keep checking my 
## ernestyalumni.wordpress.com blog and various social media channels 
## for updates as I try to keep putting out great stuff. 
## 
## Fund Science! & Help Ernest in Physics Research! : 
##  quantum super-A-polynomials - researched by Ernest Yeung
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
################################################################################
## 
#################### 
## References  
####################
# Darryl D. Holm, Tanya Schmah, Cristina Stoica. Geometric mechanics and symmetry:
#       From finite to infinite dimensions 2009

### END References
####################

SO3 = Manifold(3, 'SO3', r'SO(3)', start_index=1) ; print SO3

# SO3 is connected, but not simply connected

U = SO3.open_subset('U') ; U
cart.<x,y,z> = U.chart("x:x y:y z:z")

U_bez0    = U.open_subset('U_bez0', coord_def={cart: [x!=0,y!=0,z!=0] } )  
 # _bez0, meaning without 0 point.
cart_bez0 = cart.restrict(U_bez0)

spher.<r,ph,th> = U_bez0.chart(r'r:(0,+oo) ph:(0,2*pi):\phi th:(0,pi):\theta')

spher_to_cart = spher.transition_map( cart_bez0 , [ r * sin(th) * cos(ph) , r * sin(th) * sin(ph) , r*cos(th) ] )

spher_to_cart.set_inverse( sqrt( x^2+y^2+z^2 ) , atan2( y,x ) , atan2( sqrt( x^2 + y^2 ), z) )


cyl.<rc,phc,zc> = U_bez0.chart(r'rc:(0,+oo):r phc:(0,2*pi):\phi zc:(-oo,oo):z' )
cyl_to_cart = cyl.transition_map( cart_bez0 , ( rc*cos(phc) , rc*sin(phc) , zc ) )
cyl_to_cart.set_inverse( sqrt( x^2+y^2 ) , atan2(y,x) , z ) 


E_1 = Matrix( 3 , [ 0 ,  0 , 0 , 0 , 0 , -1 ,  0 , 1 , 0 ] )
E_2 = Matrix( 3 , [ 0 ,  0 , 1 , 0 , 0 ,  0 , -1 , 0 , 0 ] )
E_3 = Matrix( 3 , [ 0 , -1 , 0 , 1 , 0 ,  0 ,  0 , 0 , 0 ] )
Ehat_i = [ E_1, E_2, E_3 ] 

def Commutation( A, B):
    return A*B- B*A

def AntiCommutation( A,B):
    return A*B + B*A

## Test angular velocity
omega_1, omega_2, omega_3 = var("omega_1 omega_2 omega_3", domain=real)
assume(omega_1,"real")
assume(omega_2,"real")
assume(omega_3,"real")

psi_1, psi_2, psi_3 = var("psi_1 psi_2 psi_3", domain=real)
assume(psi_1,"real")
assume(psi_2,"real")
assume(psi_3,"real")

#####
## e.g. getting the rotation matrix that is simplified
sageobj( exp( omega_1 * E_1)[1][2]._maxima_().demoivre()).simplify_full()


t = var("t",domain=real)
assume(t,"real")

# test curve

omega_1t = function('omega_1t', t)
omega_2t = function('omega_2t', t)
omega_3t = function('omega_3t', t)

E_i = cart_bez0.frame()   # i = 1, 2, 3 
Ei  = cart_bez0.coframe() # i = 1, 2, 3 
omegatest  = omega_1  * E_i[1] + omega_2  * E_i[2] + omega_3  * E_i[3] 
omegatestt = omega_1t * E_i[1] + omega_2t * E_i[2] + omega_3t * E_i[3]
psitest    = psi_1  * E_i[1] + psi_2  * E_i[2] + psi_3  * E_i[3] 

def hat_map( x ):
    """
    hat_map = hat_map( x ) 
    x = ( x_1 , x_2 , x_3 )

    hat_map(x) = Matrix( [ [    0 , -x_3 ,  x_2 ] , 
    	   	     [  x_3 ,    0 , -x_1 ] ,
		     [ -x_2 ,  x_1 ,    0 ] ] )
    """
    try:
	x_1 , x_2 , x_3 = x
    	return Matrix( [ [    0 , -x_3 ,  x_2 ] , 
    	   	       	 [  x_3 ,    0 , -x_1 ] ,
		     	 [ -x_2 ,  x_1 ,    0 ] ] )
    except TypeError:
    	   x_1 , x_2, x_3 = [ x_i.expr() for x_i in x ]
    	   return Matrix( [ [    0 , -x_3 ,  x_2 ] , 
    	   	 	    [  x_3 ,    0 , -x_1 ] ,
		     	    [ -x_2 ,  x_1 ,    0 ] ] )

omegahattest = hat_map( [ omega_test[i] for i in range(1,4) ] )

Pi1, Pi2, Pi3 = var("Pi1 Pi2 Pi3", domain="real")
assume(Pi1, "real")
assume(Pi2, "real")
assume(Pi3, "real")
Pi = Pi1 * Ei[1] + Pi2 * Ei[2] + Pi3 * Ei[3] # cf. pp. 200 5.3 Isomorphisms of Lie groups and Lie algebras Example 5.36 Holm, Geometric mechanics and symmetry
Pi(omegatest).display()

print 1/2 * ( hat_map( [ Pi[i].expr() for i in range(1,4) ] ) * omegahattest.transpose()).trace()

def gen_Ad( R, xi):
    """
    gen_Ad stands for general Adjoint mapping
    gen_Ad = gen_Ad(R,xi)
    Ad: G x \mathfrak{g} -> \mathfrak{g}

    """
    return R * xi * R.inverse()

def Ad( R, omega ):
    """
    Ad = Ad( R, omega )
    It was shown on pp. 224, Holm, et. al. Geometric Mechanics, Sec. 6.2 Actions of a Lie group on itself that for the case of SO(3) and Lie Algebra so(3), the Adjoint action of SO(3) on so(3) is left multiplication action of SO(3) on R^3
    """	   
    return R*omega

def make_rotation( omega ):
    """
    make_rotation = make_rotation( omega )
    [ omega_1, omega_2, omega_3 ] = omega
    """
    Rlist = []
    for i in range(3):
    	if omega[i] == 0:
	   Rlist.append( Matrix.identity(3) )
	else:
	  R__i = exp( omega[i] * Ehat_i[i] )
	  R__icopy = []
	  for j in range(3):
	      for k in range(3):
	      	  R__icopy.append( sageobj( R__i[j][k]._maxima_().demoivre() ).simplify_full() )
 	  Rlist.append( Matrix( 3, R__icopy ) )
    return prod( Rlist )


##########
## coefficient of inertia matrix
##########

j_11 , j_12, j_13 , j_22 , j_23 , j_33 = var("j_11 j_12 j_13 j_22 j_23 j_33",domain="real")
assume(j_11,"real")
assume(j_12,"real")
assume(j_13,"real")
assume(j_22,"real")
assume(j_23,"real")
assume(j_33,"real")
J = Matrix( SR , 3 , [ j_11 , j_12, j_13, j_12 , j_22 , j_23 , j_13 , j_23 , j_33 ] )


############################################################
## Riemannian metric on SO(3), kinetic energy metric
## cf. pp. 244 Sec. 7.1 Rigid-body dynamics Holm, et. al. Geometric Mechanics
############################################################


KEm = U_bez0.riemann_metric('KEm')
KEm[:] = [ [ (hat_map( E_i[i][:] )*J*hat_map(E_i[j][:]).transpose()).trace() for j in range(1,4) ] for i in range(1,4) ] 

