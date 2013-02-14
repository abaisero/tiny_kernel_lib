// Author: Andrea Baisero
// Contact: baisero@kth.se
// Data: December 2012
//
// --------------------
// Kernel computation based on the author's own publication "The Path Kernel", ICPRAM 2012.
// --------------------

#ifndef _PATH_KERNEL_HPP_
#define _PATH_KERNEL_HPP_

#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<algorithm>
#include<vector>
#include"RefKernel.hpp"

/** @brief Path Kernel class
 *
 *  Computes the Path kernel on sequential inputs. Formally described for non-empty sequences as:
 *  \f[
 *      k_{PATH}(s,t) = k_{\Sigma}(s_1,t_1) + C_{HV} k_{PATH}(s_{2:},t) + C_{HV} k_{PATH}(s,t_{2:}) + C_D k_{PATH}(s_{2:},t_{2:})
 *  \f]
 *
 *  This kernel class extends RefKernel, and thus requires the specification of an internal kernel instance.
 *
 *  What is a sequence
 *  ------------------
 *
 *  Sequences are represented as vectors of some other data type called "symbol" type. The nature of the symbol type can be any, as long as the symbol kernel is adequate to handle it.
 *
 *  Data Inputs
 *  -----------
 *
 *  The inputs to this kernel must be vectors of whichever other type is adequate for the supplied internal instance.
 *
 *  For example, if \f$ k_{\Sigma} = k_{RBF} \f$, then the inputs to \f$ k_{PATH} \f$ must be vectors of the inputs to \f$ k_{RBF} \f$, i.e. vectors of vectors of some basic numeric type (double, float, int..).
 *
 *  How is it *really* computed?
 *  ----------------------------
 *
 *  The elaboration is based, for computational efficiency, upon an alterative expression:
 *  \f[
 *      k_{PATH}(s,t) = \sum_{i,j}^{|s||t|} k_{\Sigma}(s_i,t_i) k_{\omega}(i,j)
 *  \f]
 *  where \f$ k_{\omega} \f$ is refered to as "weight kernel" or "weight matrix" and is
 *  \f[
 *      k_{\omega}(1,1) = 1,
 *  \f]
 *  and
 *  \f[
 *      k_{\omega}(i,j) = C_{HV} k_{\omega}(i-1,j) + C_{HV} k_{\omega}(i,j-1) + C_D k_{\omega}(i-1,j-1).
 *  \f]
 *
 */
template<typename SK>
class PathKernel: public RefKernel<SK> {
	protected:
        /** @brief Cost relative to horizontal and vertical steps. */
        const double _CHV;

        /** @brief Cost relative to diagonal steps. */
        const double _CD;

        /** @brief Current weight matrix.
         *
         *  Used for efficient computations of the kernel.
         */
        std::vector<std::vector<double> > wmat;

        /** @brief Current dimension of the weight matrix. */
        size_t _DIM;

        /** @brief Directory destinated to weight-matrix files.
         *
         *  A non-empty string automatically gives read permission on the files containes in the folder.
         */
        std::string wDir;

        /** @brief Write permission in the `wDir` folder. */
        bool wW;

	public: 
        /** @brief Default value for the `_CHV` attribute. */
        static const double _CHV_def;

        /** @brief Default value for the `_CD` attribute. */
        static const double _CD_def;

        /** @brief Initializes the internal kernel reference and the step-related parameters.
         *
         *  @param[in] sk
         *          Kernel class instance. 
         *  @param[in] CHV
         *          Horizontal/Vertical step weight.
         *  @param[in] CD
         *          Diagonal step weight.
         */
        PathKernel(SK &sk,const double CHV,const double CD);

        /** @brief Initializes the internal kernel reference.
        *
        *  The step-related parameters are assigned their default values.
        *
        *  @param[in] sk
        *          Kernel class instance. 
        */
        PathKernel(SK &sk);

        /** @brief Evaluates the kernel function \f$ k_{PATH}(s,t) \f$ and stores the result in referenced parameter k.
         *
         *  @param[in] s
         *          Sequential (std::vector of symbols) input.
         *  @param[in] t 
         *          Sequential (std::vector of symbols) input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename SYM_TYPE,typename RET_TYPE>
        void operator()(const std::vector<SYM_TYPE> &s,const std::vector<SYM_TYPE> &t,RET_TYPE &k);

        /** @brief Evaluates the kernel function \f$ k_{PATH}(s,t) \f$ and stores the result in referenced parameter k.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(s,s,k)`.
         *
         *  @param[in] s
         *          Sequential (std::vector of symbols) input.
         *  @param[out] k
         *          Variable in which the kernel value is stored.
         */
        template<typename SYM_TYPE,typename RET_TYPE>
        void operator()(const std::vector<SYM_TYPE> &s,RET_TYPE &k);

        /** @brief Evaluates the kernel function \f$ k_{NORM}(s_i,t_j) \f$ with \f$ s_i\in \f$ `slist` and \f$ t_j\in \f$ `tlist`, and stores the result in reference matrix parameter km.
        *
        *  After evaluation, `km[i][j]` is set to the kernel value computed on `slist[i]` and `tlist[j]`.
        *
        *  @param[in] slist
        *          List (std::vector) of sequential (std::vector of symbols) inputs.
        *  @param[in] tlist
        *          List (std::vector) of sequential (std::vector of symbols) inputs.
        *  @param[out] km
        *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
        */
        template<typename SYM_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<SYM_TYPE> > &slist,const std::vector<std::vector<SYM_TYPE> > &tlist,std::vector<std::vector<RET_TYPE> > &km);

        /** @brief Evaluates the kernel function \f$ k_{PATH}(s_i,s_j) \f$ with \f$ s_i,s_j\in \f$ `slist`, and stores the result in reference matrix parameter km.
         *
         *  After evaluation, `km[i][j]` is set to the kernel value computed on `slist[i]` and `slist[j]`.
         *
         *  Equivalent, albeit optimised, to calling the more explicit version `(*this)(slist,slist,km)`.
         *
         *  @param[in] slist
         *          List (std::vector) of sequential (std::vector of symbols) inputs.
         *  @param[out] km
         *          Reference to a matrix (std::vector<std::vector>) variable in which the kernel values are stored.
         */
        template<typename SYM_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<SYM_TYPE> > &slist,std::vector<std::vector<RET_TYPE> > &km);

        /** @brief Evaluates the kernel function \f$ k_{PATH}(s_i,s_i) \f$ with \f$ s_i\in \f$ `slist`, and stores the result in reference vector parameter kv.
         *
         *  Equivalent to:
         *
         *      kv.resize(slist.size());
         *      for(size_t i=0;i<slist.size();i++)
         *          (*this)(slist[i],kv[i]);
         *
         *  @param[in] slist
         *          List (std::vector) of sequential (std::vector of symbols) inputs.
         *  @param[out] kv
         *          Reference to a vector (std::vector) variable in which the kernel values are stored.
         */
        template<typename SYM_TYPE,typename RET_TYPE>
        void operator()(const std::vector<std::vector<SYM_TYPE> > &slist,std::vector<RET_TYPE> &kv);

        /** @brief Updates the weight matrix to reach a specific dimension.
         *
         *  Does nothing if the weight matrix already has a dimension greater or equal to `dim`.
         *  
         *  The update operation is efficient: the current weight matrix is not re-elaborated from scratch.
         *
         *  @param[in] dim
         *          Dimension up to which to extend the weight matrix.
         */
        void updateWMat(const size_t dim);

        /** @brief Returns a copy of the weight matrix.
         *
         *  @return
         *          The weight matrix.
         */
        std::vector<std::vector<double> > getWMat();

        /** @brief Configures the kernel to enable load/save of the weight matrix.
         *
         *  @param[in] f
         *          Path to the folder where to load/save files.
         *  @param[in] w
         *          Permission to load/save files.
         *          If false, only loading is enabled. If true, saving is also enabled.
         */
        void folder(const std::string &f,const bool &w);

        /** @brief Saves current weight matrix on disk.
         *
         *  Does nothing if:
         *  - folder was not set.
         *  - write permission was not given.
         *  - an existing file, relative to the same parameters, already exists which contains a weight matrix with a dimension greater or equal to the current weight matrix.
         *
         *  @return
         *          True if a file was actually written. False otherwise (for whichever reason).
         */
        bool saveWMat() const;

        /** @brief Loads weight matrix from disk.
         *
         *  Does nothing if:
         *  - folder was not set.
         *  - a file relative to the current parameters does not exist.
         *  - the size of the current weight matrix is already greater or equal than the one in the file.
         *  
         *  @return
         *          True if the weight matrix was actually read from file. False otherwise (for whichever reason).
         */
        bool loadWMat();

    private:
        /** @brief Initializes the weight matrix to dimension 1x1.  */
        void initWMat();
}; 

template<typename SK>
const double PathKernel<SK>::_CHV_def   = 0.9/3.0;
template<typename SK>
const double PathKernel<SK>::_CD_def    = 1.1/3.0;

template<typename SK>
PathKernel<SK>::PathKernel(SK &sk,const double CHV,const double CD): RefKernel<SK>(sk),_CHV(CHV),_CD(CD),_DIM(0),wDir(""),wW(false) {
    if(CHV<=0)
        throw "Parameter \"CHV\" is not positive.";
    if(CD<=0)
        throw "Parameter \"CD\" is not positive.";
    initWMat();
};

template<typename SK>
PathKernel<SK>::PathKernel(SK &sk): RefKernel<SK>(sk),_CHV(_CHV_def),_CD(_CD_def),_DIM(0),wDir(""),wW(false) {
    initWMat();
};

template<typename SK>
void PathKernel<SK>::initWMat() {
    wmat.resize(1);
    wmat[0].resize(1);
    wmat[0][0]=1;
    _DIM=1;
}

template<typename SK>
template<typename SYM_TYPE,typename RET_TYPE>
void PathKernel<SK>::operator()(const std::vector<SYM_TYPE> &s,const std::vector<SYM_TYPE> &t,RET_TYPE &k) {
    size_t ls=s.size();
    size_t lt=t.size();
    if(ls==0||lt==0)
        return;
    updateWMat(std::max(ls,lt));
    std::vector<std::vector<RET_TYPE> > skm;
    this->_sk(s,t,skm);
    k=RET_TYPE(0);
    for(size_t i=0;i<ls;i++)
        for(size_t j=0;j<lt;j++)
            k+=skm[i][j]*(wmat[i][j]+wmat[ls-i-1][lt-j-1])/2;
}

template<typename SK>
template<typename SYM_TYPE,typename RET_TYPE>
void PathKernel<SK>::operator()(const std::vector<SYM_TYPE> &s,RET_TYPE &k) {
    size_t ls=s.size();
    if(ls==0)
        return;
    updateWMat(ls);
    std::vector<std::vector<RET_TYPE> > skm;
    this->_sk(s,skm);
    k=RET_TYPE(0);
    for(size_t i=0;i<ls;i++) {
        k+=skm[i][i]*wmat[i][i];
        for(size_t j=i+1;j<ls;j++)
            k+=2*skm[i][j]*(wmat[i][j]+wmat[ls-i-1][ls-j-1])/2;
    }
}

template<typename SK>
template<typename SYM_TYPE,typename RET_TYPE>
void PathKernel<SK>::operator()(const std::vector<std::vector<SYM_TYPE> > &slist,const std::vector<std::vector<SYM_TYPE> > &tlist,std::vector<std::vector<RET_TYPE> > &km) {
    size_t lsl=slist.size();
    size_t ltl=tlist.size();
    if(lsl==0||ltl==0)
        throw "Empty sequence vector.";
    km.resize(lsl);
    for(size_t i=0;i<lsl;i++) {
        km[i].resize(ltl);
        for(size_t j=0;j<ltl;j++)
            (*this)(slist[i],tlist[j],km[i][j]);
    }
}

template<typename SK>
template<typename SYM_TYPE,typename RET_TYPE>
void PathKernel<SK>::operator()(const std::vector<std::vector<SYM_TYPE> > &slist,std::vector<std::vector<RET_TYPE> > &km) {
    size_t lsl=slist.size();
    km.resize(lsl);
    if(lsl==0)
        throw "Empty sequence vector.";
    for(size_t i=0;i<lsl;i++) {
        km[i].resize(lsl);
        (*this)(slist[i],km[i][i]);
        for(size_t j=0;j<i;j++) {
            (*this)(slist[i],slist[j],km[i][j]);
            km[j][i]=km[i][j];
        }
    }
}

template<typename SK>
template<typename SYM_TYPE,typename RET_TYPE>
void PathKernel<SK>::operator()(const std::vector<std::vector<SYM_TYPE> > &slist,std::vector<RET_TYPE> &kv) {
    size_t lsl=slist.size();
    kv.resize(lsl);
    for(size_t i=0;i<lsl;i++)
        (*this)(slist[i],kv[i]);
}

template<typename SK>
void PathKernel<SK>::updateWMat(const size_t dim) {
    if(dim>_DIM) {
        double temp;
        size_t old_dim=_DIM;
        _DIM=dim;
        wmat.resize(_DIM);
        for(size_t i=0;i<_DIM;i++)
            wmat[i].resize(_DIM);
        for(size_t i=old_dim;i<_DIM;i++) {
            temp=_CHV*wmat[i-1][0];
            wmat[i][0]=wmat[0][i]=temp;
        }
        for(size_t i=1;i<old_dim;i++) {
            for(size_t j=old_dim;j<_DIM;j++) {
                temp=_CHV*(wmat[i-1][j]+wmat[i][j-1])+_CD*wmat[i-1][j-1];
                wmat[i][j]=wmat[j][i]=temp;
            }
        }
        for(size_t i=old_dim;i<_DIM;i++) {
            wmat[i][i]=2*_CHV*wmat[i-1][i]+_CD*wmat[i-1][i-1];
            for(size_t j=i+1;j<_DIM;j++) {
                temp=_CHV*(wmat[i-1][j]+wmat[i][j-1])+_CD*wmat[i-1][j-1];
                wmat[i][j]=wmat[j][i]=temp;
            }
        }
    }
}

template<typename SK>
std::vector<std::vector<double> > PathKernel<SK>::getWMat() {
    return wmat;
}

template<typename SK>
void PathKernel<SK>::folder(const std::string &f,const bool &w) {
    wDir=std::string(f);
    wW=bool(w);
}

template<typename SK>
bool PathKernel<SK>::saveWMat() const {
    bool save=false;
    std::string fname;
    if(wDir.size()>0 && wW) {
        std::stringstream sstr;
        sstr << wDir << "/wmat_CHV_" << std::scientific << std::setprecision(10) << _CHV << "_CD_" << _CD << ".bin";
        fname=sstr.str();
        std::ifstream ifs(fname.c_str(),std::ios::in|std::ios::binary);
        if(ifs.is_open()) {
            size_t N;
            ifs.read((char*)&N,sizeof(N));
            if(_DIM>N)
                save=true;
            ifs.close();
        }
        else
            save=true;
    }
    if(save) {
        std::ofstream ofs(fname.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
        ofs.write((char*)&_DIM,sizeof(_DIM));
        for(size_t i=0;i<_DIM;i++)
            for(size_t j=0;j<_DIM;j++)
                ofs.write((char*)&wmat[i][j],sizeof(wmat[i][j]));
        ofs.close();
    }
    return save;
}

template<typename SK>
bool PathKernel<SK>::loadWMat() {
    bool loaded=false;
    if(wDir.size()>0) {
        std::stringstream sstr;
        sstr << wDir << "/wmat_CHV_" << std::scientific << std::setprecision(10) << _CHV << "_CD_" << _CD << ".bin";
        std::string fname=sstr.str();
        std::ifstream ifs(fname.c_str(),std::ios::in|std::ios::binary);
        if(ifs.is_open()) {
            size_t N;
            ifs.read((char*)&N,sizeof(N));
            if(N>_DIM) {
                loaded=true;
                _DIM=N;
                wmat.resize(_DIM);
                for(size_t i=0;i<_DIM;i++)
                    wmat[i].resize(_DIM);
                for(size_t i=0;i<_DIM;i++)
                    for(size_t j=0;j<_DIM;j++)
                        ifs.read((char*)&wmat[i][j],sizeof(wmat[i][j]));
            }
            ifs.close();
        }
    }
    return loaded;
}

#endif // _PATH_KERNEL_HPP_

