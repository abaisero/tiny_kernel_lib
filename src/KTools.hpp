// Author: Andrea Baisero
// Contact: baisero@kth.se
// Data: December 2012
//
// --------------------
// Kernel computation based on the author's own publication "The Path Kernel", ICPRAM 2012.
// --------------------

#ifndef _KTOOLS_HPP_
#define _KTOOLS_HPP_

#include<cmath>
#include<vector>

/** @brief Namespace with useful tools for the elaboration of kernels and kernel matrices.
 *
 *  Contains methods for the elaboration of normalized kernels, normalized kernel matrices, distance matrices, and other simple matrix analysis and manipulation tools.
 *
 */
namespace ktools {

    /////////////////////////
    // NORMALIZATION TOOLS //
    /////////////////////////

    /** @brief Transforms a square kernel matrix into a normalized kernel matrix.
     *  
     *  @param[in] nkm
     *          Kernel matrix to normalize.
     */
    template<typename RET_TYPE>
    void kern2norm(std::vector<std::vector<RET_TYPE> > &nkm);

    /** @brief Evaluates the normalized kernel function \f$ \tilde k_{SK}(x,y) \f$ relative to non-normalized kernel `sk` and stores the result in reference parameter nk.
     *
     *  Alternative method to elaborate normalized kernels, which does not require to instantiate NormKernel. Defined as
     *
     *  \f$ \tilde k_{SK}(x,y) = \frac{ k_{SK}(x,y) }{ \sqrt{ k_{SK}(x,x) k_{SK}(y,y) } }\f$.
     *
     *  @param[in]  sk
     *          Kernel instance to normalize.
     *  @param[in]  x
     *          First input suitable for `sk`.
     *  @param[in]  y
     *          Second input suitable for `sk`.
     *  @param[out] nk
     *          Variable in which the normalized kernel value is stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &nk);

    /** @brief Evaluates the normalized kernel function \f$ \tilde k_{SK}(x,x) \f$ relative to non-normalized kernel `sk` and stores the result in reference parameter nk.
     *
     *  Equivalent, albeit optimised, to calling the more explicit version `norm(sk,x,x,nk)`.
     *
     *  @param[in] sk
     *          Kernel instance to normalize.
     *  @param[in] x
     *          Input suitable for `sk`.
     *  @param[out] nk
     *          Variable in which the normalized kernel value is stored.
     *          Always set to 1 unless input x is a vector of zeros.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const DATA_TYPE &x,RET_TYPE &nk);

    /** @brief Evaluates the normalized kernel function \f$ \tilde k_{SK}(x_i,y_j) \f$ with \f$ x_i\in \f$ `xlist` and \f$ y_j\in \f$ `ylist`, relative to non-normalized kernel instance `sk`, and stores the result in reference matrix parameter nkm.
     *
     *  After evaluation, `nkm[i][j]` is set to the normalized kernel value computed on `xlist[i]` and `ylist[j]`.
     *
     *  @param[in] sk
     *          Kernel instance to normalize.
     *  @param[in] xlist
     *          First list (std::vector) of inputs suitable for `sk`.
     *  @param[in] ylist
     *          Second list (std::vector) of inputs suitable for `sk`.
     *  @param[out] nkm
     *          Reference to a matrix (std::vector<std::vector>) variable in which the normalized kernel values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &nkm);

    /** @brief Evaluates the normalized kernel function \f$ \tilde k_{SK}(x_i,x_j) \f$ with \f$ x_i,x_j\in \f$ `xlist`, relative to non-normalized kernel instance `sk`, and stores the result in reference matrix parameter nkm.
     *
     *  After evaluation, `nkm[i][j]` is set to the normalized kernel value computed on `xlist[i]` and `xlist[j]`.
     *
     *  Equivalent, albeit optimised, to calling the more explicit version `ktools::norm(sk,xlist,xlist,km)`.
     *
     *  @param[in] sk
     *          Kernel instance to normalize.
     *  @param[in] xlist
     *          List (std::vector) of inputs suitable for `sk`.
     *  @param[out] nkm
     *          Reference to a matrix (std::vector<std::vector>) variable in which the normalized kernel values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &nkm);

    /** @brief Evaluates the normalized kernel \f$ \tilde k_{SK}(x_i,x_i) \f$ with \f$ x_i\in \f$ `xlist`, relative to non-normalized kernel instance `sk`, and stores the result in reference vector parameter nkv.
     *
     *  Equivalent to:
     *
     *      kv.resize(xlist.size());
     *      for(size_t i=0;i<xlist.size();i++)
     *          ktools::norm(sk,xlist[i],kv[i]);
     *
     *  @param[in] sk
     *          Kernel instance to normalize.
     *  @param[in] xlist
     *          List (std::vector) of vectorial inputs.
     *  @param[out] nkv
     *          Reference to a vector (std::vector) variable in which the normalized kernel values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &nkv);

    ////////////////////
    // DISTANCE TOOLS //
    ////////////////////

    /** @brief Transforms a square kernel matrix into a distance matrix.
     *  
     *  @param[in] dm
     *          Kernel matrix to transform.
     */
    template<typename RET_TYPE>
    void kern2dist(std::vector<std::vector<RET_TYPE> > &dm);

    /** @brief Evaluates the distance function \f$ d_{SK}(x,y) \f$ relative to kernel `sk` and stores the result in reference parameter d.
     *
     *  Distance defined as \f$ d_{SK}(x,y) = \sqrt{k_{SK}(x,x) + k_{SK}(y,y) -2k_{SK}(x,y)} \f$.
     *
     *  @param[in]  sk
     *          Kernel instance with which to elaborate the distance.
     *  @param[in]  x
     *          First input suitable for `sk`.
     *  @param[in]  y
     *          Second input suitable for `sk`.
     *  @param[out] d
     *          Variable in which the distance value is stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &d);

    /** @brief Evaluates the distance function \f$ d_{SK}(x,x) \f$ relative to kernel `sk` and stores the result in reference parameter d.
     *
     *  Equivalent, albeit optimised, to calling the more explicit version `dist(sk,x,x,d)`.
     *
     *  @param[in]  sk
     *          Kernel instance with which to elaborate the distance.
     *  @param[in]  x
     *          First input suitable for `sk`.
     *  @param[out] d
     *          Variable in which the distance value is stored.
     *          Always set to 0.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const DATA_TYPE &x,RET_TYPE &d);

    /** @brief Evaluates the distance function \f$ d_{SK}(x_i,y_j) \f$ with \f$ x_i\in \f$ `xlist` and \f$ y_j\in \f$ `ylist`, relative to kernel instance `sk`, and stores the result in reference matrix parameter dm.
     *
     *  After evaluation, `dm[i][j]` is set to the distance value computed between `xlist[i]` and `ylist[j]`.
     *
     *  @param[in] sk
     *          Kernel instance with which to elaborate the distance.
     *  @param[in] xlist
     *          First list (std::vector) of inputs suitable for `sk`.
     *  @param[in] ylist
     *          Second list (std::vector) of inputs suitable for `sk`.
     *  @param[out] dm
     *          Reference to a matrix (std::vector<std::vector>) variable in which the distance values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &dm);

    /** @brief Evaluates the distance function \f$ d_{SK}(x_i,x_j) \f$ with \f$ x_i,x_j\in \f$ `xlist`, relative to on-normalized kernel instance `sk`, and stores the result in reference matrix parameter km.
     *
     *  After evaluation, `nkm[i][j]` is set to the normalized kernel value computed on `xlist[i]` and `xlist[j]`.
     *
     *  Equivalent, albeit optimised, to calling the more explicit version `ktools::norm(sk,xlist,xlist,km)`.
     *
     *  @param[in] sk
     *          Kernel instance with which to elaborate the distance.
     *  @param[in] xlist
     *          List (std::vector) of inputs suitable for `sk`.
     *  @param[out] dm
     *          Reference to a matrix (std::vector<std::vector>) variable in which the distance values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &dm);

    /** @brief Evaluates the distance function \f$ d_{SK}(x_i,x_i) \f$ with \f$ x_i\in \f$ `xlist`, relative to kernel instance `sk`, and stores the result in reference vector parameter dv.
     *
     *  Equivalent to:
     *
     *      dv.resize(xlist.size());
     *      for(size_t i=0;i<xlist.size();i++)
     *          ktools::dist(sk,xlist[i],dv[i]);
     *
     *  @param[in] sk
     *          Kernel instance to normalize.
     *  @param[in] xlist
     *          List (std::vector) of vectorial inputs.
     *  @param[out] dv
     *          Reference to a vector (std::vector) variable in which the distance values are stored.
     */
    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &dv);

    /////////////////////////
    // OTHER GENERIC TOOLS //
    /////////////////////////

    /** @brief Resizes a matrix into \f$ NR \times NC \f$.
     *
     *  @param[in] m
     *          Matrix to resize.
     *  @param[in] NR
     *          Number of rows.
     *  @param[in] NC
     *          Number of columns.
     */
    template<typename RET_TYPE>
    void resizeMat(std::vector<std::vector<RET_TYPE> > &m,size_t NR,size_t NC);

    /** @brief Resizes a matrix to be square and into \f$ N \times N \f$.
     *
     *  Equivalent to `ktools::resizeMat(m,N,N)`.
     *  @param[in] m
     *          Matrix to resize.
     *  @param[in] N
     *          Number of rows and columns.
     */
    template<typename RET_TYPE>
    void resizeMat(std::vector<std::vector<RET_TYPE> > &m,size_t N);

    /** @brief Verifies if matrix is square.
     *
     *  @param[in] m
     *          Matrix to verify.
     *  @returns
     *          `true` if the matrix is square. `false` otherwise.
     */
    template<typename RET_TYPE>
    bool isSquare(const std::vector<std::vector<RET_TYPE> > &m);

    /** @brief Verifies if matrix is symmetric.
     *
     *  Internally also verifies that the matrix is square.
     *
     *  @param[in] m
     *          Matrix to verify.
     *  @returns
     *          `true` if the matrix is square and symmetric. `false` otherwise.
     */
    template<typename RET_TYPE>
    bool isSymmetric(const std::vector<std::vector<RET_TYPE> > &m);

    /** @brief Verifies if the distance matrix satisfies the triangle inequality.
     *
     *  The triangle inequality is defined as \f$ d(x_i,x_j) \le d(x_i,x_k) + d(x_k,x_j) \quad \forall x_i,x_j,x_k \f$.
     *
     *  Internally also verifies that the matrix is square and symmetric.
     *
     *  @param[in] dm
     *          Matrix to verify.
     *  @returns
     *          `true` if the matrix satisfies the condition. `false` otherwise.
     */
    template<typename RET_TYPE>
    bool respectsTriangleInequality(const std::vector<std::vector<RET_TYPE> > &dm);

    /** @brief Verifies if the kernel matrix satisfies the Cauchy-Schwarz inequality.
     *
     *  The Cauchy-Schwarz inequality is defined as \f$ k(x_i,x_j)^2 \le k(x_i,x_i) k(x_j,x_j) \quad \forall x_i,x_j \f$.
     *
     *  Internally also verifies that the matrix is square and symmetric.
     *
     *  @param[in] km
     *          Matrix to verify.
     *  @returns
     *          `true` if the matrix satisfies the condition. `false` otherwise.
     */
    template<typename RET_TYPE>
    bool respectsCauchySchwarz(const std::vector<std::vector<RET_TYPE> > &km);
}

namespace ktools {

    template<typename RET_TYPE>
    void kern2norm(std::vector<std::vector<RET_TYPE> > &nkm) {
        for(size_t i=0;i<nkm.size();i++) {
            for(size_t j=0;j<i;j++) {
                if(nkm[i][j]!=RET_TYPE(0)) {
                    nkm[i][j]/=std::sqrt(nkm[i][i]*nkm[j][j]);
                    nkm[j][i]=nkm[i][j];
                }
            }
        }
        for(size_t i=0;i<nkm.size();i++)
            if(nkm[i][i]!=RET_TYPE(0))
                nkm[i][i]=RET_TYPE(1);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &nk) {
        sk(x,y,nk);
        if(nk!=RET_TYPE(0)) {
            RET_TYPE xk,yk;
            sk(x,xk);
            sk(y,yk);
            nk/=std::sqrt(xk*yk);
        }
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const DATA_TYPE &x,RET_TYPE &nk) {
        sk(x,nk);
        if(nk!=RET_TYPE(0))
            nk=RET_TYPE(1);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &nkm) {
        std::vector<RET_TYPE> xkv,ykv;
        sk(xlist,ylist,nkm);
        sk(xlist,xkv);
        sk(ylist,ykv);
        for(size_t i=0;i<xkv.size();i++)
            for(size_t j=0;j<ykv.size();j++)
                if(nkm[i][j]!=RET_TYPE(0))
                    nkm[i][j]/=std::sqrt(xkv[i]*ykv[j]);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &nkm) {
        sk(xlist,nkm);
        kern2norm(nkm);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void norm(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &nkv) {
        sk(xlist,nkv);
        for(size_t i=0;i<nkv.size();i++)
            if(nkv!=RET_TYPE(0))
                nkv=RET_TYPE(1);
    }

    template<typename RET_TYPE>
    void kern2dist(std::vector<std::vector<RET_TYPE> > &dm) {
        for(size_t i=0;i<dm.size();i++)
            for(size_t j=0;j<i;j++)
                dm[i][j]=dm[j][i]=std::sqrt(dm[i][i]+dm[j][j]-dm[i][j]-dm[j][i]);
        for(size_t i=0;i<dm.size();i++)
            dm[i][i]=RET_TYPE(0);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const DATA_TYPE &x,const DATA_TYPE &y,RET_TYPE &d) {
        RET_TYPE xk,yk;
        sk(x,y,d);
        sk(x,xk);
        sk(y,yk);
        d=std::sqrt(xk+yk-2*d);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const DATA_TYPE &x,RET_TYPE &d) {
        d=RET_TYPE(0);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,const std::vector<DATA_TYPE> &ylist,std::vector<std::vector<RET_TYPE> > &dm) {
        std::vector<RET_TYPE> xkv,ykv;
        sk(xlist,ylist,dm);
        sk(xlist,xkv);
        sk(ylist,ykv);
        for(size_t i=0;i<xkv.size();i++)
            for(size_t j=0;j<ykv.size();j++)
                dm[i][j]=std::sqrt(xkv[i]+ykv[j]-2*dm[i][j]);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<std::vector<RET_TYPE> > &dm) {
        sk(xlist,dm);
        kern2dist(dm);
    }

    template<typename SK,typename DATA_TYPE,typename RET_TYPE>
    void dist(SK &sk,const std::vector<DATA_TYPE> &xlist,std::vector<RET_TYPE> &dv) {
        for(size_t i=0;i<dv.size();i++)
            dv=RET_TYPE(0);
    }

    template<typename RET_TYPE>
    void resizeMat(std::vector<std::vector<RET_TYPE> > &m,size_t NR,size_t NC) {
        m.resize(NR);
        for(size_t i=0;i<NR;i++)
            m[i].resize(NC);
    }

    template<typename RET_TYPE>
    void resizeMat(std::vector<std::vector<RET_TYPE> > &m,size_t N) {
        resizeMat(m,N,N);
    }

    template<typename RET_TYPE>
    bool isSquare(const std::vector<std::vector<RET_TYPE> > &m) {
        size_t N=m.size();
        for(size_t i=0;i<N;i++)
            if(m[i].size()!=N)
                return false;
        return true;
    }
    
    template<typename RET_TYPE>
    bool isSymmetric(const std::vector<std::vector<RET_TYPE> > &m) {
        if(!isSquare(m))
            return false;
        size_t N=m.size();
        for(size_t i=0;i<N;i++)
            for(size_t j=0;j<i;j++)
                if(m[i][j]!=m[j][i])
                    return false;
        return true;
    }

    template<typename RET_TYPE>
    bool respectsTriangleInequality(const std::vector<std::vector<RET_TYPE> > &dm) {
        if(!isSymmetric(dm))
            return false;
        size_t N=dm.size();
        for(size_t i=0;i<N;i++)
            for(size_t j=0;j<N;j++)
                for(size_t k=0;k<N;k++)
                    if(dm[i][j]>dm[i][k]+dm[k][j])
                        return false;
        return true;
    }

    template<typename RET_TYPE>
    bool respectsCauchySchwarz(const std::vector<std::vector<RET_TYPE> > &km) {
        if(!isSymmetric(km))
            return false;
        for(size_t i=0;i<km.size();i++)
            for(size_t j=0;j<i;j++)
                if(km[i][j]*km[i][j]>km[i][i]*km[j][j])
                    return false;
        return true;
    }

}

#endif // _KTOOLS_HPP_

