#ifndef _REF_KERNEL_HPP_
#define _REF_KERNEL_HPP_

/** @brief Reference Kernel class.
 *
 *  The Reference Kernel class is a base class for kernels which rely internally on some other kernel instance.
 *  It simply provides a protected attribute consisting of the internal kernel instance, and it is meant to be inherited by other classes.
 */
template<typename SK>
class RefKernel {
    protected:
        /** @brief The internal kernel instance. */
        SK &_sk;

    public:
        /** @brief Initialises the internal kernel reference.
         *
         *  @param[in] sk
         *          Kernel class instance. 
         */
        RefKernel(SK &sk);

        /** @brief Empty virtual Destructor. Good habit for base classes.
         */
        virtual ~RefKernel();

        /** @brief Reference to the internal kernel.
         *
         *  @return
         *          The reference to the internal kernel
         */
        SK& getKernelRef();
};

template<typename SK>
RefKernel<SK>::RefKernel(SK &sk): _sk(sk) {
}

template<typename SK>
RefKernel<SK>::~RefKernel() {
}

template<typename SK>
SK& RefKernel<SK>::getKernelRef() {
    return _sk;
}

#endif

