#include<cmath>
#include<iostream>
#include<vector>
#include"RbfKernel.hpp"
#include"SymKernel.hpp"
#include"PathKernel.hpp"
#include"NormKernel.hpp"
#include"KTools.hpp"
#include"Baisero.hpp"

// Only for the purpose of this test file
using std::cout;
using std::endl;
using std::fixed;
using std::left;
using std::setprecision;
using std::setw;
using std::vector;

// Base data-type for the vectorial data
typedef vector<double> InputType_Vector;
typedef size_t InputType_Index;
typedef vector<InputType_Vector> InputType_Sequence;

void usage();
void init_data();
void usage_rbf();
void usage_sym();
void usage_pathk();
void usage_normk();
void usage_kerntools();
void usage_baisero();

// for debugging purposes
void show_k();
void show_km();

// kernel value and kernel matrix used for various outputs
double k;
vector<vector<double> > km;
// vectorial data (single vectors and lists of vectors)
InputType_Vector v1,v2,v3;
vector<InputType_Vector> vlist1,vlist2;
// indexing data (single indexes and lists of indexes)
InputType_Index i1,i2,i3;
vector<InputType_Index> ilist1,ilist2;
// sequential data (single sequences and lists of sequences)
InputType_Sequence s1,s2,s3;
vector<InputType_Sequence> slist1,slist2;

int main() {
    usage();
    init_data();
    usage_rbf();
    usage_sym();
    usage_pathk();
    usage_normk();
    usage_kerntools();
    usage_baisero();

    cout << "No runtime error encountered." << endl;
    cout << "Please refer to the source code for practical examples of usage." << endl;

    return 0;
}

void usage() {
    cout << "Kernels are implemented as functors, meaning that there is a class for each type of" << endl;
    cout << "kernel which overrides operator() to elaborate the kernel values/matrices." << endl;
    cout << "The following kernels and namespaces are implemented:" << endl;
    cout << " - RbfKernel:  A gaussian/radial basis function kernel for vectorial data." << endl;
    cout << " - SymKernel:  A kernel for labeled data." << endl;
    cout << " - PathKernel: A kernel for sequential data." << endl;
    cout << " - NormKernel: A kernel which normalises other kernels." << endl;
    cout << " - ktools:     A namespace with useful functions (distance and normalisation)." << endl;
    cout << endl;
    cout << "Assume a hypothetical kernel instance, kernel inputs and lists of inputs as follows" << endl;
    cout << endl;
    cout << "\tKernelType kern;" << endl;
    cout << endl;
    cout << "\tInputType input1,input2;" << endl;
    cout << "\tinput1 = ...     // input initialization" << endl;
    cout << "\tinput2 = ...     // input initialization" << endl;
    cout << endl;
    cout << "\tvector<InputType> input_list1,input_list2;" << endl;
    cout << "\tinput_list1 = ... // input list initialization" << endl;
    cout << "\tinput_list2 = ... // input list initialization" << endl;
    cout << endl;
    cout << "Kernel outputs are returned by reference through the last parameter." << endl;
    cout << "When the inputs are single input data-types, the output is a double," << endl;
    cout << endl;
    cout << "\tdouble k;" << endl;
    cout << "\tkern(input1,input2,k);" << endl;
    cout << endl;
    cout << "When the inputs are input lists, the output is a matrix of doubles," << endl;
    cout << endl;
    cout << "\tvector<vector<double> > km;" << endl;
    cout << "\tkern(input_list1,input_list2,km);" << endl;
    cout << endl;
    cout << "When the inputs are input lists, but you are only interested in the" << endl;
    cout << "kernel evaluation of each input with itself (i.e. the diagonal of" << endl;
    cout << "km in the previous example), the output is a vector of doubles," << endl;
    cout << endl;
    cout << "\tvector<double> kv;" << endl;
    cout << "\tkern(input_list1,kv);" << endl;
    cout << endl;
    cout << "When the inputs are the same, an efficient computation of the output" << endl;;
    cout << "is available by giving one input only once," << endl;
    cout << endl;
    cout << "\tkern(input1,k);      // equivalent to kern(input1,input1,k);" << endl;
    cout << "\tkern(input_list1,km); // equivalent to kern(input_list1,input_list1,km);" << endl;
    cout << "\tkern(input_list1,kv); // equivalent to kern(input_list1,input_list1,kv);" << endl;
    cout << endl;
    cout << "NB. The type of the last parameter helps the kernel to interpret the" << endl;
    cout << "inputs as single inputs or input lists correctly." << endl;
    cout << endl;
    cout << "What follows is a practical example which is going to actually run." << endl;
    cout << endl;
}

void init_data() {
    size_t dim=3;

    cout << "\t// ================ TYPE DEFS ================= //" << endl;
    cout << endl;
    cout << "\ttypedef vector<double> InputType_Vector;" << endl;
    cout << "\ttypedef int InputType_Index;" << endl;
    cout << "\t//typedef vector<InputType_Index> InputType_Sequence; // If the symbols are interpreted through symbolic associations" << endl;
    cout << "\ttypedef vector<InputType_Vector> InputType_Sequence; // If the symbols are interpreted as vectors" << endl;
    cout << endl;

    cout << "\t// ================= INIT DATA ================= //" << endl;
    cout << endl;

    cout << "\t// vectorial data (single vectors and lists of vectors)" << endl;
    cout << "\tInputType_Vector v1,v2,v3;" << endl;
    cout << "\tvector<InputType_Vector> vlist1,vlist2;" << endl;
    cout << endl;

    v1.resize(dim);
    v2.resize(dim);
    v3.resize(dim);
    for(size_t i=0;i<dim;i++) {
        v1[i]=i;
        v2[i]=i+.5;
        v3[i]=dim-i;
    }

    vlist1.resize(2);
    vlist1[0]=v1;
    vlist1[1]=v2;
    vlist2.resize(3);
    vlist2[0]=v1;
    vlist2[1]=v2;
    vlist2[2]=v3;

    cout << "\t// indexing data (single indexes, i.e. integers, and lists of indexes)" << endl;
    cout << "\tInputType_Index i1,i2,i3;" << endl;
    cout << "\tvector<InputType_Index> ilist1,ilist2;" << endl;
    cout << endl;
    
    i1=1;
    i2=2;
    i3=3;

    ilist1.resize(2);
    ilist1[0]=i1;
    ilist1[1]=i2;
    ilist2.resize(3);
    ilist2[0]=i1;
    ilist2[1]=i2;
    ilist2[2]=i3;

    cout << "\t// sequential data (single sequences and lists of sequences)" << endl;
    cout << "\tInputType_Sequence s1,s2,s3;" << endl;
    cout << "\tvector<InputType_Sequence> slist1,slist2;" << endl;
    cout << endl;

    s1.resize(4);
    s1[0]=v1;
    s1[1]=v2;
    s1[2]=v1;
    s1[3]=v2;
    s2.resize(5);
    s2[0]=v2;
    s2[1]=v1;
    s2[2]=v2;
    s2[3]=v1;
    s2[4]=v2;
    s3.resize(6);
    s3[0]=v3;
    s3[1]=v3;
    s3[2]=v1;
    s3[3]=v2;
    s3[4]=v3;
    s3[5]=v3;

    slist1.resize(2);
    slist1[0]=s1;
    slist1[1]=s2;
    slist2.resize(3);
    slist2[0]=s1;
    slist2[1]=s2;
    slist2[2]=s3;

    cout << "\t// Data initialization is omitted, only data declaration is important." << endl;
    cout << "\t// Notice, however, that the list vectors are built upon their base instance" << endl;
    cout << "\t// in the following manner (taking sequential inputs as an example):" << endl;
    cout << "\tslist1.resize(2);" << endl;
    cout << "\tslist1[0]=s1;" << endl;
    cout << "\tslist1[1]=s2;" << endl;
    cout << "\tslist2.resize(3);" << endl;
    cout << "\tslist2[0]=s1;" << endl;
    cout << "\tslist2[1]=s2;" << endl;
    cout << "\tslist2[2]=s3;" << endl;
    cout << endl;
    cout << "\t// This justifies the fact that the kernel matrices elaborated on slist1" << endl;
    cout << "\t// result to be sub-matrices of the kernel matrices elaborated on slist1 and slist2." << endl;
    cout << endl;
} 

void usage_rbf() {
    cout << "\t// ================= USAGE RBF KERNEL ================= //" << endl;
    cout << endl;
    cout << "\t// Constructor receives standard deviation (default = 1)" << endl;
    cout << "\tRbfKernel rbf(10);" << endl;
    // Constructor receives standard deviation (default = 1)
    RbfKernel rbf(10);

    cout << "\t// Elaborates kernel value from vectorial input" << endl;
    cout << "\trbf(v1,v2,k);" << endl;
    // Elaborates kernel value from vectorial input
    rbf(v1,v2,k);
    show_k();
    cout << "\trbf(v1,k); // equal to rbf(v1,v1,k);" << endl;
    rbf(v1,k); // equal to rbf(v1,v1,k);
    show_k();

    cout << "\t// Elaborates kernel matrix values from vectors of inputs" << endl;
    cout << "\trbf(vlist1,vlist2,km);" << endl;
    // Elaborates kernel matrix values from vectors of inputs
    rbf(vlist1,vlist2,km);
    show_km();
    cout << "\trbf(vlist1,km); // equal to rbf(vlist1,vlist1,km);" << endl;
    rbf(vlist1,km); // equal to rbf(vlist1,vlist1,km);
    show_km();
}

void usage_sym() {
    cout << "\t// ================= USAGE SYM KERNEL ================= //" << endl;
    cout << endl;
    cout << "\t// Creation of custom kernel matrix" << endl;
    // Creation of custom kernel matrix
    km.resize(10);
    for(size_t i=0;i<10;i++) {
       km[i].resize(10);
        for(size_t j=0;j<10;j++)
           km[i][j]=(i+1)*(j+1);
    }
    show_km();

    cout << "\t// Constructors receives a kernel matrix" << endl;
    cout << "\tSymKernel sym1(km);" << endl;
    cout << "\t// or a matrix dimension" << endl;
    cout << "\t//SymKernel sym2(10); // equal to sym(km); where km is an Identity matrix of size 10" << endl;
    // Constructors receives a kernel matrix
    SymKernel sym1(km);
    // or a matrix dimension
    SymKernel sym2(10); // equal to sym(km); where km is an Identity matrix of size 10

    cout << "\t// Elaborates single kernel value from integer input" << endl;
    cout << "\tsym1(i1,i2,k);" << endl;
    // Elaborates single kernel value from integer input
    sym1(i1,i2,k);
    show_k();
    cout << "\tsym1(i1,k); // equal to rbf(i1,i1,k);" << endl;
    sym1(i1,k); // equal to rbf(i1,i1,k);
    show_k();

    cout << "\t// Elaborates kernel matrix values from vectors of inputs" << endl;
    cout << "\tsym1(ilist1,ilist2,km);" << endl;
    // Elaborates kernel matrix values from vectors of inputs
    sym1(ilist1,ilist2,km);
    show_km();
    cout << "\tsym1(ilist1,km); // equal to rbf(ilist1,ilist1,km);" << endl;
    sym1(ilist1,km); // equal to rbf(ilist1,ilist1,km);
    show_km();
}

void usage_pathk() {
    cout << "\t// ================= USAGE PATH KERNEL ================= //" << endl;
    cout << endl;
    cout << "\t// Creation of ground kernel" << endl;
    cout << "\tRbfKernel sk;" << endl;
    // Creation of ground kernel
    RbfKernel sk;
    cout << endl;

    cout << "\t// Constructor receives a kernel instance" << endl;
    cout << "\tPathKernel<RbfKernel> pk(sk);  // with default parameters" << endl;
    cout << "\tPathKernel<RbfKernel> pk2(sk,0.35,0.3); // with custom parameters" << endl;
    // Constructor receives a kernel instance
    PathKernel<RbfKernel> pk(sk);
    cout << endl;

    cout << "\t// Elaborates single kernel value from sequential input" << endl;
    cout << "\tpk(s1,s2,k);" << endl;
    // Elaborates single kernel value from sequential input
    pk(s1,s2,k);
    show_k();
    cout << "\tpk(s1,k); // equal to pk(s1,s1,k);" << endl;
    pk(s1,k); // equal to pk(s1,s1,k);
    show_k();

    cout << "\t// Elaborates kernel matrix values from vectors of inputs" << endl;
    cout << "\tpk(slist1,slist2,km);" << endl;
    // Elaborates kernel matrix values from vectors of inputs
    pk(slist1,slist2,km);
    show_km();
    cout << "\tpk(slist1,km); // equal to pk(slist1,slist1,km);" << endl;
    pk(slist1,km); // equal to pk(slist1,slist1,km);
    show_km();

    cout << "The PathKernel class also provides functionalities to save and load on file" << endl;
    cout << "intermediary results which may improve future performance. To activate the" << endl;
    cout << "functionality, you must provide the folder to use as storage, together with" << endl;
    cout << "a flag which determines if the class instance is allowed to write on disk" << endl;
    cout << "(by specifying a folder, read permission is automatically given)." << endl;
    cout << endl;

    cout << "\t// Increases the size of the internal matrix" << endl;
    cout << "\tint N=10;" << endl;
    cout << "\tpk.updateWMat(N); // where N is an estimate to the maximum length of the input sequences" << endl;
    cout << "\t                  // in case of doubt, there is no harm in over-estimating the value" << endl;
    cout << "\tkm=pk.getWMat();" << endl;
    // Increases the size of the internal matrix
    size_t N=10;
    pk.updateWMat(N); // where N is an estimate to the maximum length of the input sequences
                      // in case of doubt, there is no harm in over-estimating the value
    km=pk.getWMat();
    show_km();

    cout << "\t// Enables load and save on file" << endl;
    cout << "\tpk.folder(\"folder_path\",true);" << endl;
    cout << "\tpk.saveWMat();" << endl;
    cout << "\tpk.loadWMat();" << endl;
    // Enables load and save on file
    pk.folder("folder_path",true);
    pk.saveWMat();
    pk.loadWMat();
    cout << endl;
}

void usage_normk() {
    cout << "\t// ================= USAGE NORM KERNEL ================= //" << endl;
    cout << endl;
    cout << "\t// Example built upon the normalized path kernel" << endl;
    cout << "\t// Creation of path kernel" << endl;
    cout << "\tRbfKernel sk;" << endl;
    cout << "\tPathKernel<RbfKernel> pk(sk);" << endl;
    // Example built upon the normalized path kernel
    // Creation of path kernel
    RbfKernel sk;
    PathKernel<RbfKernel> pk(sk);
    cout << endl;

    cout << "\t// Constructor receives a kernel instance" << endl;
    cout << "\tNormKernel<PathKernel<RbfKernel> > nk(pk);" << endl;
    // Constructor receives a kernel instance
    NormKernel<PathKernel<RbfKernel> > nk(pk);
    cout << endl;

    cout << "\t// Elaborates single kernel value from sequential input" << endl;
    cout << "\tnk(s1,s2,k);" << endl;
    // Elaborates single kernel value from sequential input
    nk(s1,s2,k);
    show_k();
    cout << "\tnk(s1,k); // equal to nk(s1,s1,k);" << endl;
    nk(s1,k); // equal to nk(s1,s1,k);
    show_k();

    cout << "\t// Elaborates kernel matrix values from vectors of inputs" << endl;
    cout << "\tnk(slist1,slist2,km);" << endl;
    // Elaborates kernel matrix values from vectors of inputs
    nk(slist1,slist2,km);
    show_km();
    cout << "\tnk(slist1,km); // equal to nk(slist1,slist1,km);" << endl;
    nk(slist1,km); // equal to nk(slist1,slist1,km);
    show_km();
}

void usage_kerntools() {
    cout << "\t// ================= USAGE KERN TOOLS ================= //" << endl;
    cout << endl;
    cout << "The namespace \"ktools\" contains a number of functions to produce:" << endl;
    cout << " - Normalized kernel values/matrices" << endl;
    cout << " - Distance matrices" << endl;
    cout << endl;
    cout << "Both functionalities work by either receiving a kernel and the data" << endl;
    cout << "to process, or a pre-computed kernel matrix" << endl;
    cout << endl;
    
    cout << "\t// Creation of an arbitraty kernel" << endl;
    cout << "\tRbfKernel sk;" << endl;
    cout << "\tPathKernel<RbfKernel> pk(sk);" << endl;
    // Creation of an arbitraty kernel
    RbfKernel sk;
    PathKernel<RbfKernel> pk(sk);
    cout << endl;

    cout << "\t// Elaborates normalised kernel matrix" << endl;
    cout << "\tpk(slist1,km);" << endl;
    cout << "\tktools::kern2norm(km); // Only works with square kernel matrices" << endl;
    // Elaborates normalised kernel matrix
    pk(slist1,km);
    ktools::kern2norm(km); // Only works with square kernel matrices
    show_km();

    cout << "\t// Elaborates normalised kernel value" << endl;
    cout << "\tktools::norm(pk,s1,s2,k);" << endl;
    // Elaborates normalised kernel value
    ktools::norm(pk,s1,s2,k);
    show_k();
    cout << "\tktools::norm(pk,s1,k);" << endl;
    ktools::norm(pk,s1,k);
    show_k();

    cout << "\t// Elaborates normalised kernel matrix" << endl;
    cout << "\tktools::norm(pk,slist1,slist2,km);" << endl;
    // Elaborates normalised kernel matrix
    ktools::norm(pk,slist1,slist2,km);
    show_km();
    cout << "\tktools::norm(pk,slist1,km);" << endl;
    ktools::norm(pk,slist1,km);
    show_km();

    cout << "\t// Elaborates distance matrix" << endl;
    cout << "\tpk(slist1,km);" << endl;
    cout << "\tktools::kern2dist(km); // Only works with square kernel matrices" << endl;
    // Elaborates distance matrix
    pk(slist1,km);
    ktools::kern2dist(km); // Only works with square kernel matrices
    show_km();

    cout << "\t// Elaborates distance value" << endl;
    cout << "\tktools::dist(pk,s1,s2,k);" << endl;
    // Elaborates distance value
    ktools::dist(pk,s1,s2,k);
    show_k();
    cout << "\tktools::dist(pk,s1,k); // Valid, but always returns 0" << endl;
    ktools::dist(pk,s1,k); // Valid, but always returns 0
    show_k();

    cout << "\t// Elaborates normalised distance matrix" << endl;
    cout << "\tktools::dist(pk,slist1,slist2,km);" << endl;
    // Elaborates normalised distance matrix
    ktools::dist(pk,slist1,slist2,km);
    show_km();
    cout << "\tktools::dist(pk,slist1,km);" << endl;
    ktools::dist(pk,slist1,km);
    show_km();
}

void usage_baisero() {
    cout << "\t// ================= USAGE BAISERO ================= //" << endl;
    k=baisero::selectSigma(slist2,10);
    show_k();
}

void show_k() {
    cout << "\tk:" << endl;
    cout << "\t    " << left << fixed << setw(5) << setprecision(4) << k << endl << endl;
}

void show_km() {
    cout << "\tKM:" << endl;
    for(size_t i=0;i<km.size();i++) {
        cout << "\t    ";
        for(size_t j=0;j<km[i].size();j++) {
            cout << left << fixed << setw(5) << setprecision(4) << km[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

