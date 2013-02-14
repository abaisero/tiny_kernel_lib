Small Kernel Library {#mainpage}
====================

This is an *extremely* small library of C++ classes for kernel evaluations.

The main purpose of this library is to provide programmers with a computationally efficient version of the **Path Kernel**.
The kernels contained in this library are minimalistic and limited:
-   RbfKernel, computes the well-know Gaussian kernel.
-   SymKernel, computes a kernel for labeled data, based on pre-defined results.
-   NormKernel, computes normalized versions of any other kernel instance.
-   PathKernel, computes a kernel for sequential data, for a generic type of symbols.

Additionally, the following are provided:
-   RefKernel, a kernel base class for kernels which depend on other kernels (such as NormKernel and PathKernel).
-   ktools, a namespace with "useful" functions.

*Additionally*, the further 
-   baisero, a namespace with ad-hoc methods which might be required for any specific reason and are not suitable for generic applications.

What is a kernel class
----------------------

A kernel class is designed in this library to be any class which overrides the `operator()` method five times, obtaining the following method signatures:

-   `void operator()(const DATA_TYPE &x,const DATA_TYPE &y,RETURN_TYPE &k);`

    to compute the kernel value between two data-points `x` and `y`, and store the result in `k`.

-   `void operator()(const DATA_TYPE &x,RETURN_TYPE &k);`

    to compute the kernel value between a data-point `x` and itself, and store the result in `k`.

-   `void operator()(const vector<DATA_TYPE> &xlist,const vector<DATA_TYPE> &ylist,vector<vector<RETURN_TYPE> > &km);`

    to compute the kernel values between all data-points in lists `xlist` and `ylist`, and store the results in matrix `km`.

-   `void operator()(const vector<DATA_TYPE> &xlist,vector<vector<RETURN_TYPE> > &km);`

    to compute the kernel values between all pairs of data-points in list `xlist`, and store the results in matrix `km`.

-   `void operator()(const vector<DATA_TYPE> &xlist,vector<RETURN_TYPE> &kv);`

    to compute the kernel value between all data-point `x` and themselves, and store the results in vector `kv`.

DATA_TYPE may be any data type in accordance with the kernel's domain, and RETURN_TYPE should be left as a templated type.

Usage
-----

Please read the [usage](@ref usage) page.

Observations
------------

Unfortunately, because some of the proposed kernel classes and their methods are "templated", there is no base class or interface which on can use to build up kernels. It becomes thus the programmer's burden to thoroughly check that the required methods are implemented correctly.

