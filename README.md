# MATLAB Code for Mixed Finite Element Methods for Stress Gradient Elasticity

This repository contains MATLAB codes related to the numerical experiments in the paper

**Mixed Finite Element Methods for Planar Stress Gradient Elasticity**
by Ting Lin and Shudan Tian.

The codes are mainly used to reproduce the numerical experiments for the mixed finite element methods developed in the paper, including the Lagrange and Hermite-type stress finite element pairs, inf-sup tests, and convergence experiments for the stress-gradient elasticity model.

## Dependencies

The implementation is based on **iFEM**, an integrated finite element method package in MATLAB created by **Long Chen**.

The iFEM package can be found at:

```text
https://github.com/lyc102/ifem
```

Please download iFEM and add it to your MATLAB path before running the codes in this repository. For example, after downloading iFEM, you may run the following command in MATLAB from the iFEM root directory:

```matlab
setpath
```

or manually add the iFEM folder and its subfolders to the MATLAB path.

## Notes on Reproducibility

Some routines in this repository use mesh generation, finite element basis functions, quadrature rules, and degree-of-freedom management tools from iFEM. Therefore, the codes may not run correctly unless iFEM is properly installed and added to the MATLAB path.

We are currently cleaning and reorganizing the codes. In the next few days, we will also upload additional auxiliary files and selected iFEM-related routines used in our experiments, so that users who have difficulty locating or installing the original iFEM package can still reproduce the numerical results more conveniently.

## How to Run

After installing iFEM and adding it to the MATLAB path, open MATLAB and run the corresponding test scripts directly.

Typical scripts include experiments for:

* convergence tests of the proposed finite element pairs;
* inf-sup constant tests with respect to geometric parameters;
* asymptotic behavior with respect to the stress-gradient parameter;
* boundary-layer-related numerical tests.

Some scripts may require precomputed reference solutions or data files. These files will be included or documented as the repository is further cleaned.

## Acknowledgement

This code uses routines from the iFEM package developed by Long Chen. Users are encouraged to cite iFEM when using these routines.

A suggested citation for iFEM is:

```bibtex
@techreport{Chen2009iFEM,
  author = {Long Chen},
  title = {{iFEM}: An Integrated Finite Element Methods Package in {MATLAB}},
  institution = {University of California at Irvine},
  year = {2009},
  url = {https://github.com/lyc102/ifem}
}
```

## Contact

For questions about the codes or the numerical experiments, please contact shudan.tian@xtu.edu.cn
