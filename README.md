# diffgeo-by-sagemnfd
Differential Geometry by Sage Manifolds - implementations of Sage Manifolds package for Sage Math

-ernestyalumni

<h5>20150809</h5>
<b>version0.8examples</b> subdirectory has the examples from the sagemanifolds webpage 
http://sagemanifolds.obspm.fr/examples.html <br>
by verbatim.  

<h2>Rd.sage</h2>
<h3>Euclidean $\mathbb{R}^d$ Manifold in d-dimensions</h3>


```
R3 = Rd(3)
R3.equip_metric()
```

Make a time-dependent scalar field $p = p(x,y,z,t)\in \mathcal{C}^{\infty}(M)$ and a 
       time-dependent vector field $u = u(x,y,z,t) = u^i(x,y,z,t)\frac{\partial}{\partial x^i} \in \mathfrak{X}(M)$
```
p   = make_pt(R3.M)
u3t = make_ut(R3.M)
```

