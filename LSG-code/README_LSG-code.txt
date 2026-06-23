README for LSG-code
===================

MATLAB Code for Linear Stress Gradient Elasticity
-------------------------------------------------

This folder contains MATLAB programs and data files used for the numerical
experiments in the manuscript

    Mixed Finite Element Methods for Planar Stress Gradient Elasticity.

The implementation focuses on mixed finite element discretizations for the
planar linear stress gradient elasticity model. The files include routines for
degree-of-freedom numbering, local nodal basis construction, manufactured exact
solutions, right-hand side functions, boundary-layer tests, inf-sup constant
tests, and precomputed reference solution data.


Folder contents
---------------

1. Main experiment scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following files are the main scripts for the numerical experiments.

    NewestHermiteforLSGtest.m
        Main script for the Hermite-type finite element experiment.

    Newest_LagrangeforLSGtest.m
        Main script for the Lagrange-type finite element experiment.

    P7Hermite_Inf_sup_HermiteforLSGtest.m
        Script for testing the inf-sup constant for the P7 Hermite-type element.

    P7Lagrange_Inf_sup_HermiteforLSGtest.m
        Script for testing the inf-sup constant for the P7 Lagrange-type element.

    C2forLSG_iota_rate.m
        Main script for testing the convergence behavior with respect to the small parameter `iota` for the linear stress gradient elasticity problem.

2. Degree-of-freedom numbering routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files whose names start with "dof" construct local-to-global
degree-of-freedom numbering maps for different finite element spaces.

    C2dof.m
        Local-to-global degree-of-freedom numbering for the P5 H1-conforming
        element with C2 continuity at vertices.

    dofP2.m
        Local-to-global degree-of-freedom numbering for the P2 Lagrange element.

    dofP3.m
        Local-to-global degree-of-freedom numbering for the P3 Lagrange element.

    dofP4.m
        Local-to-global degree-of-freedom numbering for the P4 Lagrange element.

    dof_P4Hermite.m
        Local-to-global degree-of-freedom numbering for the P4 Hermite-type
        element.

    dof_P5Hermite.m
        Local-to-global degree-of-freedom numbering for the P5 Hermite-type
        element.

    dof_P7Hermite.m
        Local-to-global degree-of-freedom numbering for the P7 Hermite-type
        element.

    dofedgemy.m
        Auxiliary routine related to edge degrees of freedom.


3. Local nodal basis functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files whose names start with "combase_" compute local nodal basis functions.

    combase_P3Hermite.m
    combase_P4Hermite.m
    combase_P5Hermite.m
    combase_P7Hermite.m

The output of these routines gives the coefficients of the local nodal basis
functions with respect to the barycentric-coordinate polynomial basis. These
basis functions are used in element matrix assembly and in the evaluation of
finite element functions.


4. Quadrature and auxiliary routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    quadptsmy.m
        Provides quadrature points and weights.

    com_quad.m
        Auxiliary quadrature-related routine.

    Cursefunction.m
        Auxiliary function used by the numerical experiments.

    elemtau_n.m
        Auxiliary element-level routine used in the assembly or evaluation of
        stress-related quantities.


5. Manufactured solution and right-hand side for Example 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files with prefix "La_" are used for the first verification experiment.

    La_Ipr_LSG.m
        Exact or reference solution component used in Example 1.

    La_Ipr_sigma_LSG.m
        Exact or reference stress component used in Example 1.

    La_sigmax_Ipr_LSG.m
        Exact or reference derivative component related to the stress.

    La_sigmay_Ipr_LSG.m
        Exact or reference derivative component related to the stress.

    La_Ipf_LSG.m
        Right-hand side function corresponding to the manufactured exact
        solution.

Here, files with prefix "La_Ipr_" provide exact or reference solution
components, while "La_Ipf_" gives the corresponding forcing term.


6. Boundary-layer test functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files with prefix "layer_" are used for the boundary-layer experiment.

    layer_Ipr_LSG.m
        Exact or reference solution component for the boundary-layer experiment.

    layer_sigma_Ipr_LSG.m
        Exact or reference stress component for the boundary-layer experiment.

    layer_sigmax_Ipr_LSG.m
        Exact or reference derivative component related to the stress.

    layer_sigmay_Ipr_LSG.m
        Exact or reference derivative component related to the stress.

    layer_Ipf_LSG.m
        Right-hand side function for the boundary-layer test.

These functions are used to investigate the asymptotic behavior of the
numerical solution when the small parameter iota approaches zero and when
boundary-layer effects appear.


7. Precomputed reference solution data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following .mat files contain precomputed reference solution data.

    solutiota1new.mat
        Reference solution data for iota = 1.

    solutiota10e1new.mat
        Reference solution data for iota = 1.0e-01.

    solutiota10e2new.mat
        Reference solution data for iota = 1.0e-02.

    solutiota10e4new.mat
        Reference solution data for iota = 1.0e-04.

These reference solutions were computed using the P5 H1-conforming element with
C2 continuity at vertices on a mesh with size h = 1/128.

Additional data files include:

    COFla100000iota1.mat
    DGCOFla100000iota1.mat
    elemla100000iota1.mat
    nodela100000iota1.mat

These files store auxiliary coefficient, element, or node data used by the
numerical scripts.


How to run
----------

Open MATLAB and set the current working directory to this folder, for example

    cd LSG-code
    addpath(genpath(pwd))

Example 1: verification experiments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following two scripts correspond to the verification experiments.

    NewestHermiteforLSGtest

and

    Newest_LagrangeforLSGtest

They can be run directly in MATLAB. Before running, please make sure that the
value of the small parameter iota in the script is consistent with the
precomputed reference solution data loaded in the script.

The matching between iota and the reference data file is as follows.

    iota = 1          <-->    solutiota1new.mat
    iota = 1.0e-01    <-->    solutiota10e1new.mat
    iota = 1.0e-02    <-->    solutiota10e2new.mat
    iota = 1.0e-04    <-->    solutiota10e4new.mat

For example, if the script is run with iota = 1.0e-02, then the corresponding
reference data should be solutiota10e2new.mat. If iota is changed in the script,
the loaded reference data file should be changed accordingly.

Example 2: inf-sup constant tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following scripts are used for the numerical tests of the inf-sup constant.

    P7Hermite_Inf_sup_HermiteforLSGtest

and

    P7Lagrange_Inf_sup_HermiteforLSGtest

These scripts test the behavior of the discrete inf-sup constant for the P7
Hermite-type and P7 Lagrange-type elements, respectively.

Example 3: iota rate tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Iota-rate test

The script

```
C2forLSG_iota_rate.m
```

is used to test the convergence behavior with respect to the small parameter `iota` for the linear stress gradient elasticity problem.

This script uses the P5 H1-conforming stress element with C2 continuity at vertices and a fixed mesh size

```
hsize = 1/64.
```

The value of the small parameter is set inside the script, for example,

```
iota = 1.0e-03;
```

To test different values of `iota`, modify this line and run the script again.

#### How to run

Open MATLAB, set the current working directory to `LSG-code`, and run

```
C2forLSG_iota_rate
```

directly.

The script does not require input arguments. After the computation, the result will be printed in the MATLAB command window. The output has the form

```
Hsize:= ...
iota = ..., energy norm error = ...
```

The printed energy norm error is used to study the convergence rate with respect to `iota`.

#### Choice of exact solution and right-hand side

There are two sets of manufactured solutions in the code.

1. Boundary-layer case

This case corresponds to the test where the limiting solution does not satisfy the additional boundary condition. It is used to observe the boundary-layer effect.

Use the following functions:

```
layer_Ipf_LSG.m
    Right-hand side.

layer_Ipr_LSG.m
    Exact displacement.

layer_sigma_Ipr_LSG.m
    Exact stress.

layer_sigmax_Ipr_LSG.m
    x-derivative of the exact stress.

layer_sigmay_Ipr_LSG.m
    y-derivative of the exact stress.
```

In the current version of `C2forLSG_iota_rate.m`, this boundary-layer case is selected.

2. Non-boundary-layer case

This case corresponds to the compatible test where the limiting solution satisfies the additional boundary condition. It is used to observe the higher-order asymptotic behavior with respect to `iota`.

Use the following functions:

```
La_Ipf_LSG.m
    Right-hand side.

La_Ipr_LSG.m
    Exact displacement.

La_Ipr_sigma_LSG.m
    Exact stress.

La_sigmax_Ipr_LSG.m
    x-derivative of the exact stress.

La_sigmay_Ipr_LSG.m
    y-derivative of the exact stress.
```

To switch between the two cases, the right-hand side and all exact-solution evaluations in the error computation should be changed consistently. For example, when using the non-boundary-layer case, replace the `layer_*` functions by the corresponding `La_*` functions.


Notes
-----

1. The code is research-oriented and is intended to reproduce the numerical
   experiments reported in the manuscript.

2. File names are kept close to the notation and implementation used in the
   experiments.

3. The "dof*" routines define local-to-global numbering maps for the
   corresponding finite element spaces.

4. The "combase_*" routines compute local nodal basis functions and output their
   coefficients in the barycentric-coordinate polynomial basis.

5. The "La_*" files correspond to Example 1, namely the manufactured-solution
   verification test.

6. The "layer_*" files correspond to the boundary-layer experiment.

7. The "solutiota*new.mat" files contain reference numerical solutions for
   different values of the small parameter iota.
