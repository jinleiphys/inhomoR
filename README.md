# inhomoR
**InhomoR** is a computer code develop by Jin Lei and Pierre Descouvemont to solve the inhomogeneous equations. We applied three different methods, namely **R-matrix method**, **Green's function method**, and **Numerov method**. 

## Input description
### &global namelist:  Ecm,hcm,rmax,lmin,lmax,nr
- Ecm : Center of mass energy 
- hcm : radial step when using Numerov method.
- rmax: maximum number of the radial distance
- lmin: minimum number of partial wave 
- lmax: maximum number of partial wave
- nr:grid number when using R-matrix method and Green's function method 


### &systems: massa, massb, za,zb,ja,jb
- mass,massb : mass numbers of interaction particles 
- za, zb: charge numbers 
- ja,jb: spins 

### &potential namelist:  ptype, a1, a2, rc, uv, av, rv, uw, aw, rw, vsov, rsov, asov, vsow, rsow, asow, vd, avd, rvd, wd, awd, rwd
- a1, a2, rc, uv, av, rv, uw, aw, rw, vsov, rsov, asov, vsow, rsow, asow, vd, avd, rvd, wd, awd, rwd : potential parameters
- ptype: defines different types of potential

      1. Woods-Saxon potential.  Potential parameters are defined through the variables: 
         * Volume WS: uv, av, rv, uw, aw, rw
         * Spin-orbit: vsov, rsov, asov, vsow, rsow, asow
         * Surface (derivative) WS: vd, avd, rvd, wd, awd, rwd
      2. Gaussian potential, defined through the parameters: uv,rv,av,uw,rw,aw
