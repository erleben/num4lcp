This library contains a toolbox of Matlab functions that implements numerical methods for computing solutions for linear complementarity problems (LCPs).

The intend with this library is that users can learn and study the algorithms behind the numerical methods before selecting their favorite method and porting it to whatever target language or platform they wish for. The hope is that this will allow other researchers and practitioners to easily adopt new numerical methods into their own projects.

| Method | A-matrix properties | Available implementation | GPU capable | Subsolvers | Globalization strategies |
|:-------|:--------------------|:-------------------------|:------------|:-----------|:-------------------------|
| Minimum map Newton | No assumptions`*` | Matlab, C++/CUSP, Python | Yes (C++/CUSP) | _N/A_ | Projected Armijo backtracking line search |
| Fischer-Burmeister Newton | No assumptions | Matlab, Python C++/CUSP | Yes (C++/CUSP) | Random, Zero, Perturbation, Approximation, Penalized Chen, Chen, and Kanzow style`**` | Projected Armijo backtracking line search |
| Projected Gauss-Seidel/Projected Successive Over-Relaxation (QP) | Symmetric PSD | Matlab, C++/CUSP | No | _N/A_ | _N/A_ |
| Interior Point | No assumptions | Matlab, Python | No | Central Trajectory`***` | _N/A_ |
| Multilevel | Symmetric PD | Matlab | No | _N/A_ | _N/A_ |
| Quadratic Program | Symmetric | Matlab | No | _N/A_ | _N/A_ |

`*` C++/CUSP implementation requires A-matrix to be PSD, because it uses PCG to find search direction.

`**` Matlab implements all five subsolvers. Python implements: Random, Perturbation and Zero. C++/CUSP implements: Perturbation

`***` Newton subsystem is current solved with a Moore Penrose pseudo-inverse if badly conditioned and a direct solver otherwise (Matlab's in-built)

The context for the numerical methods is physics-based animation which means the library has only been verified and tested on typical LCPs encountered in physics-based animation.

The library was developed for academic and educational reasons. Therefore Matlab was chosen. If one does not have Matlab accessible then one may be able to easily convert the Matlab code into NumPhy, SciPhy or Octave.

Please use this reference when citing this project
```
@misc{num4lcp.11,
  author = 	 {Erleben, K. and Andersen, M. and Niebe S. and Silcowitz M.},
  title = 	 {num4lcp},
  howpublished = {Published online at code.google.com/p/num4lcp/},
  month = 	 {October},
  year = 	 2011,
  note = 	 {Open source project for numerical methods for
                  linear complementarity problems in physics-based animation}
}
```

Documentation can be found in these SIGGRAPH course notes http://iphys.wordpress.com/2013/05/19/numerical-methods-for-linear-complementarity-problems-in-physics-based-animation/