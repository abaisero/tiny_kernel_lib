#ifndef _RBF_KERNEL_HPP_
#define _RBF_KERNEL_HPP_

#include<cmath>
#include<vector>

/** @brief Radial Basis Function Kernel class
 *
 *  Computes the radial basis kernel, a.k.a Gaussian kernel, on vectorial inputs. Formally described as:
 *  \f[
 *      k_{RBF}(x,y) = e^{-\frac{\|x-y\|^2}{2\sigma^2}}
 *  \f]
 *
 *  Data Inputs
 *  -----------
 *
 *  The inputs to this kernel must be STD vectors of some basic numeric type (double, float, int..).
 *
 *  The input vectors must have the same length within the same call of the kernel functions.
 *  The length may however change in between kernel evaluation calls.
 */
class RbfKernel {

    protected:
        /** @brief Parameter relative to the radial basis function's variance.
         *
         *  Initialized to tsigma = \f$ -\frac{1}{2\sigma^2} \f$, where \f$\sigma\f$ is the user-defined standard deviation.
         */
        const double tsigma;

    public:
        /** @brief Initiates tsigma to \f$ -\frac{1}{2\sigma^2} \f$.
         *
         *  @param[in] sigma
         *          Standard deviation of the rbf (defaults to 1).
         */
        RbfKernel(double sigma=1.0);

        /** @brief Empty virtual Destructor. Good habit for base classes. */
        virtual ~RbfKernel();

        /** @brief Evaluates the kernel function \f$ k_{RBF}(x,y) \f$ and stores the result in reference parameter k.
         *
         *  @param[in] x
         *          Vectorial input.
         *  @param[in] y 
         *          Vectorial input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename VEC_TYPE,typename RET_TYPE>
        void operator()(const std::vector<VEC_TYPE> &x,const std::vector<VEC_TYPE> &y,RET_TYPE &k) const;

        /** @brief Evaluates the kernel function \f$ k_{RBF}(x,x) \f$ and stores the result in reference parameter k.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(x,x,k)`.
         *
         *  @param[in] x
         *          Vectorial input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         *          Always set to 1 unless input x is a vector of zeros.
         */
        template<typename VEC_TYPE,typename RET_TYPE>
        void operator()(const std::vector<VEC_TYPE> &x,RET_TYPE &k) const;

        /** @brief Evaluates the kernel function \f$ k_{RBF}(x_i,y_j) \f$ with \f$ x_i\in \f$ `xlist` and \f$ y_j\in \f$ `ylist`, and stores the result in reference matrix parameter km.
         *
         *  After evaluation, `km[i][j]` is set to the kernel value computed on `xlist[i]` and `ylist[j]`.
         *
         *  @param[in] xlist
         *          List (std::vector) of vectorial inputs.
         *  @param[in] ylist
         *          List (std::vector) of vectorial inputs.
         *  @param[out] km
         *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
         */
        template<typename VEC_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,const std::vector<std::vector<VEC_TYPE> > &ylist,std::vector<std::vector<RET_TYPE> > &km) const;

        /** @brief Evaluates the kernel function \f$ k_{RBF}(x_i,x_j) \f$ with \f$ x_i,x_j\in \f$ `xlist`, and stores the result in reference matrix parameter km.
         *
         *  After evaluation, `km[i][j]` is set to the kernel value computed on `xlist[i]` and `xlist[j]`.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(xlist,xlist,km)`.
         *
         *  @param[in] xlist
         *          List (std::vector) of vectorial inputs.
         *  @param[out] km
         *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
         */
        template<typename VEC_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,std::vector<std::vector<RET_TYPE> > &km) const;

        /** @brief Evaluates the kernel function \f$ k_{RBF}(x_i,x_i) \f$ with \f$ x_i\in \f$ `xlist`, and stores the result in reference vector parameter kv.
         *
         *  Equivalent to:
         *
         *      kv.resize(xlist.size());
         *      for(size_t i=0;i<xlist.size();i++)
         *          (*this)(xlist[i],kv[i]);
         *
         *  @param[in] xlist
         *          List (std::vector) of vectorial inputs.
         *  @param[out] kv
         *          Reference to a vector (std::vector) variable in which the kernel values are stored.
         */
        template<typename VEC_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,std::vector<RET_TYPE> &kv) const;
};

RbfKernel::RbfKernel(double sigma): tsigma(-1/(2*sigma*sigma)) {
    if(sigma==0)
        throw "Input \"sigma\" is 0.";
}

RbfKernel::~RbfKernel() {
}

template<typename VEC_TYPE,typename RET_TYPE>
void RbfKernel::operator()(const std::vector<VEC_TYPE> &x,const std::vector<VEC_TYPE> &y,RET_TYPE &k) const {
    if(x.empty()||y.empty())
        throw "Input vector is empty.";
    if(x.size()!=y.size())
        throw "Input vectors do not have equal size.";
    double sq_norm=0;
    for(size_t i=0;i<x.size();i++)
        sq_norm+=(x[i]-y[i])*(x[i]-y[i]);
    k=RET_TYPE(exp(tsigma*sq_norm));
}

template<typename VEC_TYPE,typename RET_TYPE>
void RbfKernel::operator()(const std::vector<VEC_TYPE> &x,RET_TYPE &k) const {
    if(x.empty())
        throw "Input vector is empty.";
    k=RET_TYPE(0);
    for(size_t i=0;i<x.size();i++) {
        if(x[i]!=RET_TYPE(0)) {
            k=RET_TYPE(1);
            break;
        }
    }
}

template<typename VEC_TYPE,typename RET_TYPE>
void RbfKernel::operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,const std::vector<std::vector<VEC_TYPE> > &ylist,std::vector<std::vector<RET_TYPE> > &km) const {
    size_t lxl=xlist.size();
    size_t lyl=ylist.size();
    if(lxl==0||lyl==0)
        throw "Input set doesn't contain any vector.";
    size_t dim=xlist[0].size();
    for(size_t i=0;i<lxl;i++) {
        if(xlist[i].empty())
            throw "Input vector is empty.";
        if(xlist[i].size()!=dim)
            throw "Input vectors do not have equal size.";
    }
    km.resize(lxl);
    for(size_t i=0;i<lxl;i++) {
        km[i].resize(lyl);
        for(size_t j=0;j<lyl;j++)
            (*this)(xlist[i],ylist[j],km[i][j]);
    }
}

template<typename VEC_TYPE,typename RET_TYPE>
void RbfKernel::operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,std::vector<std::vector<RET_TYPE> > &km) const {
    size_t lxl=xlist.size();
    if(lxl==0)
        throw "Input set doesn't contain any vector.";
    size_t dim=xlist[0].size();
    for(size_t i=0;i<lxl;i++) {
        if(xlist[i].size()==0)
            throw "Input vector is empty.";
        if(xlist[i].size()!=dim)
            throw "Input vectors do not have equal size.";
    }
    km.resize(lxl);
    for(size_t i=0;i<lxl;i++) {
        km[i].resize(lxl);
        (*this)(xlist[i],km[i][i]);
        for(size_t j=0;j<i;j++) {
            (*this)(xlist[i],xlist[j],km[i][j]);
            km[j][i]=km[i][j];
        }
    }
}

template<typename VEC_TYPE,typename RET_TYPE>
void RbfKernel::operator()(const std::vector<std::vector<VEC_TYPE> > &xlist,std::vector<RET_TYPE> &kv) const {
    size_t lxl=xlist.size();
    if(lxl!=0)
        throw "Input set doesn't contain any vector.";
    size_t dim=xlist[0].size();
    for(size_t i=0;i<lxl;i++) {
        if(xlist[i].size()==0)
            throw "Input vector is empty.";
        if(xlist[i].size()!=dim)
            throw "Input vectors do not have equal size.";
    }
    kv.resize(lxl);
    for(size_t i=0;i<lxl;i++)
        (*this)(xlist[i],kv[i]);
}

#endif // _RBF_KERNEL_HPP_

