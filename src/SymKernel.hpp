#ifndef _SYM_KERNEL_HPP_
#define _SYM_KERNEL_HPP_

#include<vector>

/** @brief Symbolic Kernel class
 *
 *  The Symbolic Kernel class provides kernel evaluation for indexing/labeled data. Formally described as:
 *  \f[
 *      k_{SYM}(i,j) = SKM_{ij}
 *  \f]
 *  where \f$ SKM \f$ is an externally-provided positive semi-definite matrix.
 *
 *  Data Inputs
 *  -----------
 *
 *  The inputs to this kernel must be integers in between 0 (included) and the dimension of \f$ SKM \f$ (excluded).
 */
class SymKernel {

    protected:
        /** @brief Kernel matrix which characterizes the indexes/labels' kernel values.
         *  
         *  Although the characteristic kernel matrix may belong to any basic data type, the kernel values are internally stored as doubles.
         *  A conversion to double is executed at construction. A further conversion back to a provided basic type is executed on kernel evaluation.
         */
        std::vector<std::vector<double> > _skm;
        /** @brief Dimension of the characteristic kernel matrix. */
        size_t _N;

    public:
        /** @brief Initiates the characteristic kernel matrix
         *
         *  @param[in] skm
         *          The characteristic kernel matrix to use.
         */
        template<typename VAL_TYPE>
        SymKernel(const std::vector<std::vector<VAL_TYPE> > &skm);

        /** @brief Initiates the characteristic kernel matrix to be the identity matrix.
         *
         *  @param[in] N
         *          Dimension of the characteristic kernel matrix.
         */
        SymKernel(const size_t N);

        /** @brief Empty virtual Destructor. Good habit for base classes. */
        virtual ~SymKernel() {};

        /** @brief Evaluates the kernel function \f$ k_{SYM}(ii,jj) \f$ and stores the result in reference parameter k.
         *
         *  @param[in] ii
         *          Indexing/labeled input.
         *  @param[in] jj 
         *          Indexing/labeled input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename RET_TYPE>
        void operator()(const size_t ii,const size_t jj,RET_TYPE &k) const;

        /** @brief Evaluates the kernel function \f$ k_{SYM}(ii,ii) \f$ and stores the result in reference parameter k.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(ii,ii,k)`.
         *
         *  @param[in] ii
         *          Indexing/labeled input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename RET_TYPE>
        void operator()(const size_t ii,RET_TYPE &k) const;

        /** @brief Evaluates the kernel function \f$ k_{SYM}(ii,jj) \f$ with \f$ ii\in \f$ `ilist` and \f$ jj\in \f$ `jlist`, and stores the result in reference matrix parameter km.
         *
         *  After evaluation, `km[i][j]` is set to the kernel value computed on `ilist[i]` and `jlist[j]`.
         *
         *  @param[in] ilist
         *          List (std::vector) of indexing/labeled input.
         *  @param[in] jlist
         *          List (std::vector) of indexing/labeled input.
         *  @param[out] km
         *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
         */
        template<typename RET_TYPE>
        void operator()(const std::vector<size_t> &ilist,const std::vector<size_t> &jlist,std::vector<std::vector<RET_TYPE> > &km) const;

        /** @brief Evaluates the kernel function \f$ k_{SYM}(ii,jj) \f$ with \f$ ii,jj\in \f$ `ilist`, and stores the result in reference matrix parameter km.
         *
         *  After evaluation, `km[i][j]` is set to the kernel value computed on `ilist[i]` and `ilist[j]`.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(ilist,ilist,km)`.
         *
         *  @param[in] ilist
         *          List (std::vector) of indexing/labeled input.
         *  @param[out] km
         *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
         */
        template<typename RET_TYPE>
        void operator()(const std::vector<size_t> &ilist,std::vector<std::vector<RET_TYPE> > &km) const;

        /** @brief Evaluates the kernel function \f$ k_{SYM}(ii,ii) \f$ with \f$ ii\in \f$ `ilist`, and stores the result in reference vector parameter kv.
         *
         *  Equivalent to:
         *
         *      kv.resize(ilist.size());
         *      for(size_t i=0;i<ilist.size();i++)
         *          (*this)(ilist[i],kv[i]);
         *
         *  @param[in] ilist
         *          List (std::vector) of indexing/labeled inputs.
         *  @param[out] kv
         *          Reference to a vector (std::vector) variable in which the kernel values are stored.
         */
        template<typename RET_TYPE>
        void operator()(const std::vector<size_t> &ilist,std::vector<RET_TYPE> &kv) const;
};

template<typename VAL_TYPE>
SymKernel::SymKernel(const std::vector<std::vector<VAL_TYPE> > &skm): _skm(skm) {
    _N=skm.size();
    if(_N==0)
        throw "Parameter \"skm\" is empty.";
    for(size_t i=0;i<_N;i++) {
        if(skm[i].size()!=_N)
            throw "Parameter \"skm\" is not a square matrix.";
        for(size_t j=0;j<i;j++)
            if(skm[i][j]!=skm[j][i])
                throw "Parameter \"skm\" is not a symmetric matrix.";
    }
    _skm.resize(_N);
    for(size_t i=0;i<_N;i++) {
        _skm[i].resize(_N);
        for(size_t j=0;j<i;j++)
            _skm[i][j]=_skm[j][i]=double(skm[i][j]);
    }
}

SymKernel::SymKernel(const size_t N): _N(N) {
    if(N==0)
        throw "Parameter \"N\" is not positive.";
    _skm.resize(N);
    for(size_t i=0;i<N;i++) {
        _skm[i].resize(N,0);
        _skm[i][i]=1;
    }
}

template<typename RET_TYPE>
void SymKernel::operator()(const size_t ii,const size_t jj,RET_TYPE &k) const {
    if(ii<0||jj<0)
        throw "Input kernel index is negative.";
    if(ii>=_N||jj>=_N)
        throw "Input kernel index exceeds maximum value.";
    k=RET_TYPE(_skm[ii][jj]);
}

template<typename RET_TYPE>
void SymKernel::operator()(const size_t ii,RET_TYPE &k) const {
    if(ii<0)
        throw "Input kernel index is negative.";
    if(ii>=_N)
        throw "Input kernel index exceeds maximum value.";
    k=RET_TYPE(_skm[ii][ii]);
}

template<typename RET_TYPE>
void SymKernel::operator()(const std::vector<size_t> &ilist,const std::vector<size_t> &jlist,std::vector<std::vector<RET_TYPE> > &km) const {
    size_t lil=ilist.size();
    size_t ljl=jlist.size();
    if(lil==0||ljl==0)
        throw "Empty kernel index vector.";
    for(size_t i=0;i<lil;i++) {
        if(ilist[i]<0)
            throw "Input kernel index is negative.";
        if(ilist[i]>=_N)
            throw "Input kernel index exceeds maximum value.";
    }
    for(size_t j=0;j<ljl;j++) {
        if(jlist[j]<0)
            throw "Input kernel index is negative.";
        if(jlist[j]>=_N)
            throw "Input kernel index exceeds maximum value.";
    }
    km.resize(lil);
    for(size_t i=0;i<lil;i++) {
        km[i].resize(ljl);
        for(size_t j=0;j<ljl;j++)
            km[i][j]=RET_TYPE(_skm[ilist[i]][jlist[j]]);
    }
}

template<typename RET_TYPE>
void SymKernel::operator()(const std::vector<size_t> &ilist,std::vector<std::vector<RET_TYPE> > &km) const {
    size_t lil=ilist.size();
    if(lil==0)
        throw "Empty kernel index vector.";
    for(size_t i=0;i<lil;i++) {
        if(ilist[i]<0)
            throw "Input kernel index is negative.";
        if(ilist[i]>=_N)
            throw "Input kernel index exceeds maximum value.";
    }
    km.resize(lil);
    for(size_t i=0;i<lil;i++) {
        km[i].resize(lil);
        km[i][i]=RET_TYPE(_skm[ilist[i]][ilist[i]]);
        for(size_t j=0;j<i;j++)
            km[i][j]=km[j][i]=RET_TYPE(_skm[ilist[i]][ilist[j]]);
    }
}

template<typename RET_TYPE>
void SymKernel::operator()(const std::vector<size_t> &ilist,std::vector<RET_TYPE> &kv) const {
    size_t lil=ilist.size();
    if(lil==0)
        throw "Empty kernel index vector.";
    for(size_t i=0;i<lil;i++) {
        if(ilist[i]<0)
            throw "Input kernel index is negative.";
        if(ilist[i]>=_N)
            throw "Input kernel index exceeds maximum value.";
    }
    kv.resize(lil);
    for(size_t i=0;i<lil;i++)
        kv[i]=RET_TYPE(_skm[ilist[i]][ilist[i]]);
}

#endif // _SYM_KERNEL_HPP_

