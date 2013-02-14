// Author: Andrea Baisero
// Contact: baisero@kth.se
// Data: December 2012
//
// --------------------
// Kernel computation based on the author's own publication "The Path Kernel", ICPRAM 2012.
// --------------------

#ifndef _NORM_KERNEL_HPP_
#define _NORM_KERNEL_HPP_

#include<cmath>
#include<vector>

/** @brief Normalized Kernel class.
 *
 *  Computes the normalized version of some other kernel \f$ SK \f$. Formally described as:
 *  \f[
 *      k_{NORM}(x,y) = \frac{k_{SK}(x,y)}{\sqrt{k_{SK}(x,x)k_{SK}(y,y)}}
 *  \f]
 *
 *  This kernel class extends RefKernel, and thus requires the specification of an internal kernel instance.
 *
 *  Data Inputs
 *  -----------
 *
 *  The inputs to this kernel must be adequate for the supplied internal kernel instance.
 */
template<typename SK>
class NormKernel: public RefKernel<SK>{
	public: 
        /** @brief Initiates the internal kernel.
         *
         *  @param[in] sk
         *          Kernel to normalize.
         */
        NormKernel(SK &sk);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(x,y) \f$ and stores the result in reference parameter k.
         *
         *  @param[in] x
         *          Vectorial input.
         *  @param[in] y 
         *          Vectorial input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename DATA_TYPE,typename RET_TYPE>
        void operator()(const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &k);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(x,x) \f$ and stores the result in reference parameter k.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(x,x,k)`.
         *
         *  @param[in] x
         *          Vectorial input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         *          Always set to 1 unless input x is a vector of zeros.
         */
        template<typename DATA_TYPE,typename RET_TYPE>
        void operator()(const DATA_TYPE &x,RET_TYPE &k);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(x_i,y_j) \f$ with \f$ x_i\in \f$ `xlist` and \f$ y_j\in \f$ `ylist`, and stores the result in reference matrix parameter km.
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
        template<typename DATA_TYPE,typename RET_TYPE>
        void operator()(const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &km);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(x_i,x_j) \f$ with \f$ x_i,x_j\in \f$ `xlist`, and stores the result in reference matrix parameter km.
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
        template<typename DATA_TYPE,typename RET_TYPE>
        void operator()(const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &km);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(x_i,x_i) \f$ with \f$ x_i\in \f$ `xlist`, and stores the result in reference vector parameter kv.
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
        template<typename DATA_TYPE,typename RET_TYPE>
        void operator()(const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &kv);
}; 

template<typename SK>
NormKernel<SK>::NormKernel(SK &sk): RefKernel<SK>(sk) {
};

template<typename SK>
template<typename DATA_TYPE,typename RET_TYPE>
void NormKernel<SK>::operator()(const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &k) {
    RET_TYPE kx,ky;
    this->_sk(x,y,k);
    if(k!=RET_TYPE(0)) {
        this->_sk(x,x,kx);
        this->_sk(y,y,ky);
        k=RET_TYPE(k/std::sqrt(kx*ky));
    }
}

template<typename SK>
template<typename DATA_TYPE,typename RET_TYPE>
void NormKernel<SK>::operator()(const DATA_TYPE &x,RET_TYPE &k) {
    this->_sk(x,k);
    if(k!=RET_TYPE(0))
        k=RET_TYPE(1);
}

template<typename SK>
template<typename DATA_TYPE,typename RET_TYPE>
void NormKernel<SK>::operator()(const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &km) {
    std::vector<RET_TYPE> kvx,kvy;
    this->_sk(xlist,ylist,km);
    this->_sk(xlist,kvx);
    this->_sk(ylist,kvy);
    for(size_t i=0;i<xlist.size();i++)
        for(size_t j=0;j<ylist.size();j++)
            if(km[i][j]!=RET_TYPE(0))
                km[i][j]=RET_TYPE(km[i][j]/std::sqrt(kvx[i]*kvy[j]));
}

template<typename SK>
template<typename DATA_TYPE,typename RET_TYPE>
void NormKernel<SK>::operator()(const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &km) {
    this->_sk(xlist,km);
    for(size_t i=0;i<xlist.size();i++)
        for(size_t j=0;j<i;j++)
            if(km[i][j]!=RET_TYPE(0))
                km[i][j]=km[j][i]=RET_TYPE(km[i][j]/std::sqrt(km[i][i]*km[j][j]));
    for(size_t i=0;i<xlist.size();i++)
        if(km[i][i]!=RET_TYPE(0))
            km[i][i]=RET_TYPE(1);
}

template<typename SK>
template<typename DATA_TYPE,typename RET_TYPE>
void NormKernel<SK>::operator()(const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &kv) {
    this->_sk(xlist,kv);
    for(size_t i=0;i<xlist.size();i++)
        if(kv[i]!=RET_TYPE(0))
            kv[i]=RET_TYPE(1);
}

#endif // _NORM_KERNEL_HPP_

