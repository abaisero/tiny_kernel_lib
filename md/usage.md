Usage Example {#usage}
=============

This page explains the basics at using some of the kernels in this library.
To run all the examples, check out the [code](@ref usage_code).

Generic Usage rules
---------------------

Assume a hypothetical kernel instance, kernel inputs and lists of inputs as follows

~~~{.c}
KernelType kern;

InputType input1,input2;
input1 = ...     // input initialization
input2 = ...     // input initialization

vector<InputType> input_list1,input_list2;
input_list1 = ... // input list initialization
input_list2 = ... // input list initialization
~~~

Kernel outputs are returned by reference through the last parameter.
When the inputs are single input data-types, the output is a double,

~~~{.c}
double k;
kern(input1,input2,k);
~~~

When the inputs are input lists, the output is a matrix of doubles,

~~~{.c}
vector<vector<double> > km;
kern(input_list1,input_list2,km);
~~~

When the inputs are input lists, but you are only interested in the
kernel evaluation of each input with itself (i.e. the diagonal of
km in the previous example), the output is a vector of doubles,

~~~{.c}
vector<double> kv;
kern(input_list1,kv);
~~~

When the inputs are the same, an efficient computation of the output"
is available by giving one input only once,

~~~{.c}
kern(input1,k);      // equivalent to kern(input1,input1,k);
kern(input_list1,km); // equivalent to kern(input_list1,input_list1,km);
kern(input_list1,kv); // equivalent to kern(input_list1,input_list1,kv);
~~~

NB. The type of the last parameter helps the kernel to interpret the
inputs as single inputs or input lists correctly.

What follows is a practical example which is going to actually run.

Typedefs
--------

~~~{.c}
typedef vector<double> InputType_Vector;
typedef int InputType_Index;
//typedef vector<InputType_Index> InputType_Sequence; // If the symbols are interpreted through symbolic associations
typedef vector<InputType_Vector> InputType_Sequence; // If the symbols are interpreted as vectors
// Sequences can actually be vectors of *ANY* other data-type.
~~~

Data Types and Initialization
-----------------------------

~~~{.c}
// vectorial data (single vectors and lists of vectors)
InputType_Vector v1,v2,v3;
vector<InputType_Vector> vlist1,vlist2;

// indexing data (single indexes, i.e. integers, and lists of indexes)
InputType_Index i1,i2,i3;
vector<InputType_Index> ilist1,ilist2;

// sequential data (single sequences and lists of sequences)
InputType_Sequence s1,s2,s3;
vector<InputType_Sequence> slist1,slist2;

// Data initialization is omitted, only data declaration is important.
// Notice, however, that the list vectors are built upon their base instance
// in the following manner (taking sequential inputs as an example):
slist1.resize(2);
slist1[0]=s1;
slist1[1]=s2;
slist2.resize(3);
slist2[0]=s1;
slist2[1]=s2;
slist2[2]=s3;

// This justifies the fact that the kernel matrices elaborated on slist1
// result to be sub-matrices of the kernel matrices elaborated on slist1 and slist2.

~~~

Usage RbfKernel
---------------

~~~{.c}
// Constructor receives standard deviation (default = 1)
RbfKernel rbf(10);
// Elaborates kernel value from vectorial input
rbf(v1,v2,k);
k:
    0.9963

rbf(v1,k); // equal to rbf(v1,v1,k);
k:
    1.0000

// Elaborates kernel matrix values from vectors of inputs
rbf(vlist1,vlist2,km);
KM:
    1.0000 0.9963 0.9465 
    0.9963 1.0000 0.9572 

rbf(vlist1,km); // equal to rbf(vlist1,vlist1,km);
KM:
    1.0000 0.9963 
    0.9963 1.0000 
~~~

Usage SymKernel
---------------

~~~{.c}
// Creation of custom kernel matrix
KM:
    1.0000 2.0000 3.0000 4.0000 5.0000 6.0000 7.0000 8.0000 9.0000 10.0000 
    2.0000 4.0000 6.0000 8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000 
    3.0000 6.0000 9.0000 12.0000 15.0000 18.0000 21.0000 24.0000 27.0000 30.0000 
    4.0000 8.0000 12.0000 16.0000 20.0000 24.0000 28.0000 32.0000 36.0000 40.0000 
    5.0000 10.0000 15.0000 20.0000 25.0000 30.0000 35.0000 40.0000 45.0000 50.0000 
    6.0000 12.0000 18.0000 24.0000 30.0000 36.0000 42.0000 48.0000 54.0000 60.0000 
    7.0000 14.0000 21.0000 28.0000 35.0000 42.0000 49.0000 56.0000 63.0000 70.0000 
    8.0000 16.0000 24.0000 32.0000 40.0000 48.0000 56.0000 64.0000 72.0000 80.0000 
    9.0000 18.0000 27.0000 36.0000 45.0000 54.0000 63.0000 72.0000 81.0000 90.0000 
    10.0000 20.0000 30.0000 40.0000 50.0000 60.0000 70.0000 80.0000 90.0000 100.0000 


// Constructors receives a kernel matrix
SymKernel sym1(km);
// or a matrix dimension
SymKernel sym2(10); // equal to sym(km); where km is an Identity matrix of size 10
// Elaborates single kernel value from integer input
sym1(i1,i2,k);
k:
    6.0000

sym1(i1,k); // equal to rbf(i1,i1,k);
k:
    4.0000

// Elaborates kernel matrix values from vectors of inputs
sym1(ilist1,ilist2,km);
KM:
    4.0000  6.0000  8.0000  
    6.0000  9.0000  12.0000

sym1(ilist1,km); // equal to rbf(ilist1,ilist1,km);
KM:
    4.0000  6.0000  
    6.0000  9.0000  
~~~

Usage PathKernel
----------------

~~~{.c}
// Creation of ground kernel
RbfKernel sk;

// Constructor receives a kernel instance
PathKernel<RbfKernel> pk(sk);  // with default parameters
PathKernel<RbfKernel> pk2(sk,0.35,0.3); // with custom parameters

// Elaborates single kernel value from sequential input
pk(s1,s2,k);
k:
    4.1211

pk(s1,k); // equal to pk(s1,s1,k);
k:
    3.8949

// Elaborates kernel matrix values from vectors of inputs
pk(slist1,slist2,km);
KM:
    3.8949 4.1211 1.5171 
    4.1211 4.8300 1.8324 

pk(slist1,km); // equal to pk(slist1,slist1,km);
KM:
    3.8949 4.1211 
    4.1211 4.8300 
~~~

The PathKernel class also provides functionalities to save and load on file
intermediary results which may improve future performance. To activate the
functionality, you must provide the folder to use as storage, together with
a flag which determines if the class instance is allowed to write on disk
(by specifying a folder, read permission is automatically given).

~~~{.c}
// Increases the size of the internal matrix
int N=10;
pk.updateWMat(N); // where N is an estimate to the maximum length of the input sequences
                  // in case of doubt, there is no harm in over-estimating the value
km=pk.getWMat();
KM:
    1.0000 0.3000 0.0900 0.0270 0.0081 0.0024 0.0007 0.0002 0.0001 0.0000 
    0.3000 0.5467 0.3010 0.1314 0.0518 0.0192 0.0069 0.0024 0.0008 0.0003 
    0.0900 0.3010 0.3810 0.2641 0.1429 0.0676 0.0294 0.0121 0.0047 0.0018 
    0.0270 0.1314 0.2641 0.2982 0.2292 0.1414 0.0760 0.0372 0.0170 0.0074 
    0.0081 0.0518 0.1429 0.2292 0.2468 0.2005 0.1348 0.0795 0.0426 0.0212 
    0.0024 0.0192 0.0676 0.1414 0.2005 0.2108 0.1772 0.1265 0.0799 0.0459 
    0.0007 0.0069 0.0294 0.0760 0.1348 0.1772 0.1836 0.1580 0.1177 0.0784 
    0.0002 0.0024 0.0121 0.0372 0.0795 0.1265 0.1580 0.1621 0.1419 0.1092 
    0.0001 0.0008 0.0047 0.0170 0.0426 0.0799 0.1177 0.1419 0.1446 0.1282 
    0.0000 0.0003 0.0018 0.0074 0.0212 0.0459 0.0784 0.1092 0.1282 0.1299 

// Enables load and save on file
pk.folder("folder_path",true);
pk.saveWMat();
pk.loadWMat();
~~~

Usage NormKernel
----------------

~~~{.c}
// Example built upon the normalized path kernel
// Creation of path kernel
RbfKernel sk;
PathKernel<RbfKernel> pk(sk);

// Constructor receives a kernel instance
NormKernel<PathKernel<RbfKernel> > nk(pk);

// Elaborates single kernel value from sequential input
nk(s1,s2,k);
k:
    0.9501

nk(s1,k); // equal to nk(s1,s1,k);
k:
    1.0000

// Elaborates kernel matrix values from vectors of inputs
nk(slist1,slist2,km);
KM:
    1.0000 0.9501 0.3738 
    0.9501 1.0000 0.4054 

nk(slist1,km); // equal to nk(slist1,slist1,km);
KM:
    1.0000 0.9501 
    0.9501 1.0000 
~~~

Usage KTools
------------

The namespace "KTools" contains a number of functions to produce:
 - Normalized kernel values/matrices
 - Distance matrices

Both functionalities work by either receiving a kernel and the data
to process, or a pre-computed kernel matrix

~~~{.c}
// Creation of an arbitraty kernel
RbfKernel sk;
PathKernel<RbfKernel> pk(sk);

// Elaborates normalised kernel matrix
pk(slist1,km);
KTools::kern2norm(km); // Only works with square kernel matrices
KM:
    1.0000 0.9501 
    0.9501 1.0000 

// Elaborates normalised kernel value
KTools::norm(pk,s1,s2,k);
k:
    0.9501

KTools::norm(pk,s1,k);
k:
    1.0000

// Elaborates normalised kernel matrix
KTools::norm(pk,slist1,slist2,km);
KM:
    1.0000 0.9501 0.3738 
    0.9501 1.0000 0.4054 

KTools::norm(pk,slist1,km);
KM:
    1.0000 0.9501 
    0.9501 1.0000 

// Elaborates distance matrix
pk(slist1,km);
KTools::kern2dist(km); // Only works with square kernel matrices
KM:
    0.0000 0.6948 
    0.6948 0.0000 

// Elaborates distance value
KTools::dist(pk,s1,s2,k);
k:
    0.6948

KTools::dist(pk,s1,k); // Valid, but always returns 0
k:
    0.0000

// Elaborates normalised distance matrix
KTools::dist(pk,slist1,slist2,km);
KM:
    0.0000 0.6948 2.2562 
    0.6948 0.0000 2.3226 

KTools::dist(pk,slist1,km);
KM:
    0.0000 0.6948 
    0.6948 0.0000 
~~~

