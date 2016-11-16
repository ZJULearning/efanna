//#define char16_t LIBRARY_char16_t
#include <mex.h>
//#undef char16_t

template<class base> class handle_wrapper
{
public:
    handle_wrapper(base *ptr) : handle(ptr) {}
    ~handle_wrapper() {delete handle;}
    base *get_handle() {return handle;}
private:
    base *handle;
};

template<class base> inline mxArray *handle2mat(base *ptr)
{
    mxArray *mat = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(mat)) = reinterpret_cast<uint64_t>(new handle_wrapper<base>(ptr));
    return mat;
}

template<class base> inline handle_wrapper<base> *mat2wrapper(const mxArray *mat)
{
    if (mxGetNumberOfElements(mat) != 1 || mxGetClassID(mat) != mxUINT64_CLASS || mxIsComplex(mat))
        mexErrMsgTxt("Input must be a real uint64 scalar.");
    handle_wrapper<base> *wrapper = reinterpret_cast<handle_wrapper<base> *>(*((uint64_t *)mxGetData(mat)));
    return wrapper;
}

template<class base> inline base *mat2handle(const mxArray *mat)
{
    return mat2wrapper<base>(mat)->get_handle();
}

template<class base> inline void delete_wrapper(const mxArray *mat)
{
    delete mat2wrapper<base>(mat);
}
