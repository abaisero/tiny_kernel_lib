#ifndef _BAISERO_HPP_
#define _BAISERO_HPP_

#include<cmath>
#include<cstdlib>
#include<vector>

/** @brief Namespace with ad-hoc and non-general-purpose tools, which may be useful only in specific settings.
 */
namespace baisero {

    /** @brief Selects Sigma parameter for a RbfKernel when used inside a PathKernel.
     *
     *  This method receives a list of sequences (which are going to be ideally fed into the PathKernel)
     *  and outputs an educated estimate of which \f$ \sigma \f$ to use for the underlying RbfKernel.
     *
     *  Purpose:
     *  -   Useful only when using the PathKernel on the RbfKernel.
     *  -   Created for Marianna.
     *
     *  @param[in] slist
     *          A list of sequences.
     *  @param[in] N
     *          The number of iterations for the educated guess.
     *          Defaults to 0.
     *          If non-positive, it is derived from the data itself.
     *  @returns
     *          The proposed \f$\sigma\f$.
     */
    template<typename SYM_TYPE>
    double selectSigma(const std::vector<std::vector<std::vector<SYM_TYPE> > > &slist,int N=0);
}

namespace baisero {
    template<typename SYM_TYPE>
    double selectSigma(const std::vector<std::vector<std::vector<SYM_TYPE> > > &slist,int N) {
        std::vector<double> symDist(N);
        int sA,sB,symA,symB;
        size_t listSize=slist.size();
        size_t symDim=slist[0][0].size();

        if(N<=0) {
            N=0;
            for(size_t i=0;i<slist.size();i++) {
                N+=slist[i].size();
            }
            N=std::floor(std::sqrt(N));
        }

        for(int i=0;i<N;i++) {
            sA=rand()%listSize;
            symA=rand()%slist[sA].size();

            sB=rand()%listSize;
            symB=rand()%slist[sA].size();
            
            SYM_TYPE sum=SYM_TYPE(0);
            for(size_t d=0;d<symDim;d++) {
                sum+=(slist[sA][symA][d]-slist[sB][symB][d])*(slist[sA][symA][d]-slist[sB][symB][d]);
            }
            symDist[i]=std::sqrt(sum);
        }
        std::vector<double>::iterator first=symDist.begin();
        std::vector<double>::iterator last=symDist.end();
        std::vector<double>::iterator median=first+(last-first)/2;
        nth_element(first,median,last);
        return *median;
    }

}

#endif // _BAISERO_HPP_
