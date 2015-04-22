## Schwarzchild_BH.sage
## This is my implementation of Schwarzschild spacetime from 
## Eric Gourgoulhon, Michal Bejger (2014) 
##   using these authors' sagemanifolds 0.7 package for Sage Math
## The main reference that I'll liberally copy from is the pdf 
## SM_Schwarzschild on the authors' webpage 
##   for sagemanifolds
## cf. Lecture 22: Black Holes (International Winter School on Gravity and Light 2015)
############################################################################ 
## Copyleft 2015, Ernest Yeung <ernestyalumni@gmail.com>                            
##                                                            
## 20150421
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
## If you like what I'm doing and would like to help and contribute support,     
## please take a look at my crowdfunding campaign at ernestyalumni.tilt.com 
## and Patreon  
## read my mission statement and give your financial support, no matter how 
## small or large, if you can and to keep checking my 
## ernestyalumni.wordpress.com blog and various social media channels                 
## for updates as I try to keep putting out great stuff.                     
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
############################################################################ 
## 20150421 EY :  I am implementing  
## Lecture 22 Black Holes (Schwarzschild)
##   WE Heraeus International Winter School on Gravity and Light
##

M = Manifold(4, 'M', r'M' ) ; print M # default start index is 0 which is what you want

#####
## Physical constants (i.e. parameters)
#####
G_N = var('G_N') ; assume(G_N >= 0)
M_0 = var('M_0') ; assume(M_0 >= 0)
m = var('m') ; assume(m>0)

###############
## Charts on M
###############

regI  = M.open_subset('R_I', r'\mathcal{R}_{\mathrm{I}}')  # asymptotically flat regions outside the event horizon

regII = M.open_subset('R_II', r'\mathcal{R}_{\mathrm{II}}') # inside the future event horizon

regIII = M.open_subset('R_III', r'\mathcal{R}_{\mathrm{III}}') # outside event horizon
regIV = M.open_subset('R_IV', r'\mathcal{R}_{\mathrm{IV}}')

##########
#### Boyer-Lindquist coordinates
## Schwarzschild coordinates, or standard Boyer-Lindquist coordinates
##########

regI_II = regI.union(regII) ; regI_II
X.<t,r,th,ph> = regI_II.chart(r't r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')

X_I = X.restrict(regI, r>2*m)

X_II = X.restrict(regII, r<2*m)


###############
## Metric Structures
###############

g = M.lorentz_metric('g')

g[0,0] , g[1,1] = -(1-2*m/r), 1/(1-2*m/r)
g[2,2] , g[3,3] = r^2, (r*sin(th))^2


#####
## Levi-Civita connection
#####
nab = g.connection(); print nab

#####
## Curvature
#####

R = g.riemann()

#####
## Ricci tensor
#####

Ric = g.ricci()

########################################
## Eddington-Finkelstein coordinates
########################################

X_EF.<v,r,th,ph> = regI_II.chart(r'v r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\varphi')

ch_BL_EF_I = X_I.transition_map(X_EF, [t+r+2*m*ln( r/(2*m) - 1 ), r , th , ph ] , restrictions2=r>2*m )

X_EF_I = X_EF.restrict(regI) ; X_EF_I

ch_BL_EF_II = X_II.transition_map(X_EF , [t+r+2*m*ln(1-r/(2*m)), r , th, ph ] , restrictions2=r<2*m )

X_EF_II = X_EF.restrict(regII)

ch_EF_BL_I  = ch_BL_EF_I.inverse() 
ch_EF_BL_II = ch_BL_EF_II.inverse()

gI = regI.lorentz_metric('gI')

gI[0,0] , gI[1,1] = -(1-2*m/r), 1/(1-2*m/r)
gI[2,2] , gI[3,3] = r^2, (r*sin(th))^2

gII = regII.lorentz_metric('gII')

gII[0,0] , gII[1,1] = -(1-2*m/r), 1/(1-2*m/r)
gII[2,2] , gII[3,3] = r^2, (r*sin(th))^2

#####
## radial null rays
## cf. http://en.wikipedia.org/wiki/Eddington–Finkelstein_coordinates
#####

X_EF_null.<tbar,r,th,ph> = regI_II.chart(r'tbar r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\varphi')

ch_BL_EF_I_null = X_I.transition_map(X_EF_null, [t+2*m*ln( r/(2*m) - 1 ), r , th , ph ] , restrictions2=r>2*m )

X_EF_I_null = X_EF_null.restrict(regI) ; X_EF_I

ch_BL_EF_II_null = X_II.transition_map(X_EF_null , [t+2*m*ln(1-r/(2*m)), r , th, ph ] , restrictions2=r<2*m )

X_EF_II_null = X_EF_null.restrict(regII)

ch_EF_BL_I_null  = ch_BL_EF_I_null.inverse() 
ch_EF_BL_II_null = ch_BL_EF_II_null.inverse()


################################################################################ 
## Plot of the Boyer-Lindquist coordinates in terms of the Eddington-Finkelstein ones
################################################################################

X_I.plot(X_EF_I, ranges={t:(0,8) , r:(2.1,10)}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(r,v), style={t:'--', r:'-'}, parameters={m:1})

plot1 = X_I.plot(X_EF_I, ranges={t:(0,8) , r:(2.1,10)}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(r,v), style={t:'--', r:'-'}, parameters={m:1})

# plot1.show()

########################################
## Kruskal-Szekeres coordinates
########################################

M0 = regI.union(regII).union(regIII).union(regIV) ; M0

X_KS.<U,V,th,ph> = M0.chart(r'U V th:(0,pi):\theta ph:(0,2*pi):\varphi')
X_KS.add_restrictions( V^2 < 1 + U^2 )

X_KS_I = X_KS.restrict( regI, [U>0, V<U, V>-U]) ; X_KS_I

ch_BL_KS_I = X_I.transition_map(X_KS_I, [sqrt(r/(2*m)-1)*exp(r/(4*m))*cosh(t/(4*m)), sqrt((r/2*m)-1)*exp(r/(4*m))*sinh(t/(4*m)), th, ph])

X_KS_II = X_KS.restrict(regII, [V>0,V>U,V>-U]) ; X_KS_II

ch_BL_KS_II = X_II.transition_map( X_KS_II, [ sqrt(1-r/(2*m))*exp(r/(4*m))*sinh(t/(4*m)), sqrt(1-r/(2*m))*exp(r/(4*m))*cosh(t/(4*m)), th, ph] )

################################################################################
## Plot of the Boyer-Lindquist coordinates in terms of the Kruskal ones
################################################################################

graphI = X_I.plot(X_KS, ranges={t:(-12,12) , r:(2.001,5)}, nb_values={t:17, r:9}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(U,V), style={t:'--',r:'-'}, parameters={m:1})
graphII = X_II.plot(X_KS, ranges={t:(-12,12), r:(0.001,1.999)}, nb_values={t:17, r:9}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(U,V), style={t:'--', r:'-'}, color='green', parameters={m:1})
# show(graphI+graphII, xmin=-3,xmax=3, ymin=-3,ymax=3,axes_labels=['$U$','$V$'])

######################################################################
## 20150421 EY :  I am implementing concepts presented in  
## Lecture 22 Black Holes
##   WE Heraeus International Winter School on Gravity and Light
######################################################################

g.display()

################################################################################
## Curves in $M$ (or open region $U$) implemented as a class object in Python 
################################################################################

class SM_Curve(object):
    """
    SM_Curve - Sage Manifold Curve
    Initialize as such:
    SM_Curve = SM_Curve( M, U, chart )
     where M is manifold, U is open subset

    e.g. (Example of Usage)
    vee_gamma = SM_Curve( M, X.domain(), X )


    """
    def __init__(self, mnfd, U, chart ):
        self.mnfd  = mnfd
        self.U     = U
        self.chart = chart 

        self.tau = var("tau")
        assume(self.tau,"real")

        self.gamma_cur = [ function('gamma' + str(mu), self.tau) for mu in range( self.mnfd.dim())]
        self.dotgamma = [self.gamma_cur[mu].diff(self.tau) for mu in range( self.mnfd.dim())]
        self.dotgamma = sum( [ self.dotgamma[mu] * self.chart.frame()[mu] for mu in range( self.mnfd.dim() ) ] )

        self.v = [ var('v'+str(mu),domain="real") for mu in range( self.mnfd.dim() ) ]
        for mu in range( self.mnfd.dim() ):
            assume( self.v[mu], "real")
        self.v_gamma = sum( [ self.v[mu]*self.chart.frame()[mu] for mu in range(self.mnfd.dim() ) ] )        

######################################################################
## Example of use of SM_Curve and
## around minute 12 of Radial null geodesics of Lecture 22 Black Holes
## International Winter School on Gravity and Light 2015
######################################################################

# Let
nullgeo = SM_Curve( M, X.domain() , X )

# The integrand inside the action S is simply

g( nullgeo.v_gamma , nullgeo.v_gamma ).expr()

# Looks unclear, try this:
g( nullgeo.v_gamma, nullgeo.v_gamma).expr().coefficient(v0^2).full_simplify()
g( nullgeo.v_gamma, nullgeo.v_gamma).expr().coefficient(v1^2).full_simplify()
g( nullgeo.v_gamma, nullgeo.v_gamma).expr().coefficient(v2^2).full_simplify()
g( nullgeo.v_gamma, nullgeo.v_gamma).expr().coefficient(v3^2).full_simplify()

# Then the expression for the integrand of action S was indeed obtained (for "East Coast" convention of the metric)

# t-eqn. of motion
g( nullgeo.v_gamma, nullgeo.v_gamma ).expr().diff(v0).full_simplify()

## r-eqn. of motion for radial null geodesics, so theta, phi constant
## g( nullgeo.v_gamma, nullgeo.v_gamma ).expr(),diff(v1).full_simplify()

k = var('k',domain="real")
assume(k,"real")

# \dot{t}
solve( g( nullgeo.v_gamma, nullgeo.v_gamma ).expr().diff(v0).full_simplify() == k,v0)[0].rhs()

integrate( solve( g( nullgeo.v_gamma, nullgeo.v_gamma ).expr().diff(v0).full_simplify() == k,v0)[0].rhs(), r) # -1/2*(2*m*log(-2*m + r) + r)*k # so modulo -k/2 constant, the \widetilde{t}(r) is reobtained


pic1_01 = parametric_plot( [ r , ( r + 2*m * ln( r - 2 * m)).subs(m==1) ], (r, 2*1+0.25, 3.5), rgbcolor=hue(0.2) ) 

pic1_02 = parametric_plot( [ r , ( r + 2*m * ln( abs(r - 2 * m)) ).subs(m==1) ], (r, 0, 2*1-0.25), rgbcolor=hue(0.25) ) 

pic1_03 = parametric_plot( [ r , -( r + 2*m * ln( r - 2 * m)).subs(m==1) ], (r, 2*1+0.25, 3.5), rgbcolor=hue(0.95) ) 

pic1_04 = parametric_plot( [ r , -( r + 2*m * ln( abs(r - 2 * m)) ).subs(m==1) ], (r, 0, 2*1-0.25), rgbcolor=hue(0.8) ) 


show( pic1_01+pic1_02 +pic1_03+pic1_04)


## Eddington-Finkelstein
#pic2_I = X_I.plot(X_EF_I, ranges={t:(0,0.5) , r:(2.1,10)}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(r,v), style={t:'--', r:'-'}, parameters={m:1}, steps={t:0.1 , r:0.25})

pic2_I = X_I.plot(X_EF_I, ranges={tbar:(0,8) , r:(2.1,10)}, fixed_coords={th:pi/2,ph:0}, ambient_coords=(r,tbar), style={t:'--', r:'-'}, parameters={m:1})


gI.display( X_EF_I_null.frame() )
gII.display( X_EF_II_null.frame() )



