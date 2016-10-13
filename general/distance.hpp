#ifndef EFANNA_DISTANCE_H_
#define EFANNA_DISTANCE_H_
#include <stdlib.h>
#include <vector>
#include <set>
#include <xmmintrin.h>
#include <cmath>
#include <iostream> //for debug

#ifdef __GNUC__
#ifdef __AVX__
#define KGRAPH_MATRIX_ALIGN 32
#else
#ifdef __SSE2__
#define KGRAPH_MATRIX_ALIGN 16
#else
#define KGRAPH_MATRIX_ALIGN 4
#endif
#endif
#endif



namespace efanna {

template<typename T>
struct Candidate {
    size_t row_id;
    T distance;
    Candidate(const size_t row_id, const T distance): row_id(row_id), distance(distance) { }

    bool operator >(const Candidate& rhs) const {
        if (this->distance == rhs.distance) {
            return this->row_id > rhs.row_id;
        }
        return this->distance > rhs.distance;
    }
    bool operator <(const Candidate& rhs) const {
        if (this->distance == rhs.distance) {
            return this->row_id < rhs.row_id;
        }
        return this->distance < rhs.distance;
    }
};

template<typename T>
class Distance {
public:
    virtual float compare(const T* a, const T* b, size_t length) const = 0;
    virtual ~Distance() {}
};

#define SSE_L2SQR(addr1, addr2, dest, tmp1, tmp2) \
    tmp1 = _mm_load_ps(addr1);\
    tmp2 = _mm_load_ps(addr2);\
    tmp1 = _mm_sub_ps(tmp1, tmp2); \
    tmp1 = _mm_mul_ps(tmp1, tmp1); \
    dest = _mm_add_ps(dest, tmp1);


template<typename T>
class L2DistanceSSE: public Distance<T> {
public:
    typedef T ResultType;
    /**
     * Copied from flann
     * We use msse intrinstic here, we should ensure data align.
     */
ResultType compare(const T* a, const T* b, size_t size) const {
//float float_l2sqr_sse2 (float const *t1, float const *t2, unsigned dim) {
    __m128 sum;
    __m128 l0, l1, l2, l3;
    __m128 r0, r1, r2, r3;
    unsigned D = (size + 3) & ~3U;
    unsigned DR = D % 16;
    unsigned DD = D - DR;
    const float *l = a;
    const float *r = b;
    const float *e_l = l + DD;
    const float *e_r = r + DD;
    float unpack[4] __attribute__ ((aligned (16))) = {0, 0, 0, 0};
    ResultType ret = 0.0;
    sum = _mm_load_ps(unpack);
    switch (DR) {
        case 12:
            SSE_L2SQR(e_l+8, e_r+8, sum, l2, r2);
        case 8:
            SSE_L2SQR(e_l+4, e_r+4, sum, l1, r1);
        case 4:
            SSE_L2SQR(e_l, e_r, sum, l0, r0);
    }
    for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
        SSE_L2SQR(l, r, sum, l0, r0);
        SSE_L2SQR(l + 4, r + 4, sum, l1, r1);
        SSE_L2SQR(l + 8, r + 8, sum, l2, r2);
        SSE_L2SQR(l + 12, r + 12, sum, l3, r3);
    }
    _mm_storeu_ps(unpack, sum);
    ret = unpack[0] + unpack[1] + unpack[2] + unpack[3];
    return ret;//sqrt(ret);
}


    /*
    ResultType compare(const T* a, const T* b, size_t size) const {
        ResultType result = ResultType();
        size_t blocklength = size / 16 * 16; //cycle unroll 4
        __m128 sum;sum = _mm_setzero_ps();
        if(size>16)sseL2Dist(a,b,sum,blocklength);//std::cout<<blocklength<<"bl"<<std::endl;
        if(blocklength != size){
	    size_t remainder = blocklength;
	    for(;remainder < size; remainder += 4){
              __m128 tmpa = _mm_load_ps(a+remainder);
              __m128 tmpb = _mm_load_ps(b+remainder);
	      __m128 tmp = _mm_sub_ps(tmpa, tmpb);
	      tmp = _mm_mul_ps(tmp, tmp);
	      sum = _mm_add_ps(sum, tmp);
             
	    }//no more remaider cause aligned with 4 * 32 bits
        }
        ResultType* res = (ResultType*)&sum;
        for(int i = 0; i < 4; i++){
          result += res[i];
        }
        
        return result;//sqrt(result);
    }

*/
    void sseL2Dist(const T* a, const T* b, __m128& vsum, size_t length) const {
      ResultType result = ResultType();
      __m128* va = (__m128*)a;
      __m128* vb = (__m128*)b;
      //std::cout<<length<<std::endl;
      for(size_t i = 0; i < length/16*4; i += 4){
        __m128 tmp = _mm_sub_ps(va[i], vb[i]);
        tmp = _mm_mul_ps(tmp, tmp);
        vsum = _mm_add_ps(vsum, tmp);
        tmp = _mm_sub_ps(va[i+1], vb[i+1]);
        tmp = _mm_mul_ps(tmp, tmp);
        vsum = _mm_add_ps(vsum, tmp);
        tmp = _mm_sub_ps(va[i+2], vb[i+2]);
        tmp = _mm_mul_ps(tmp, tmp);
        vsum = _mm_add_ps(vsum, tmp);
        tmp = _mm_sub_ps(va[i+3], vb[i+3]);
        tmp = _mm_mul_ps(tmp, tmp);
        vsum = _mm_add_ps(vsum, tmp);
      }
     
    }
};

template<typename T>
class L2Distance: public Distance<T> {
public:
    typedef T ResultType;
    /**
     * Copied from flann
     * We do not want msse intrinstic here to avoid misalign problems.
     */
    ResultType compare(const T* a, const T* b, size_t size) const {
        ResultType result = ResultType();
        ResultType diff0, diff1, diff2, diff3;
        const T* last = a + size;
        const T* lastgroup = last - 3;

        /* Process 4 items with each loop for efficiency. */
        while (a < lastgroup) {
            diff0 = (ResultType)(a[0] - b[0]);
            diff1 = (ResultType)(a[1] - b[1]);
            diff2 = (ResultType)(a[2] - b[2]);
            diff3 = (ResultType)(a[3] - b[3]);
            result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
            a += 4;
            b += 4;
        }
        /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
        while (a < last) {
            diff0 = (ResultType)(*a++ - *b++);
            result += diff0 * diff0;
        }
        return result;//sqrt(result);
    }
};

}
#endif
