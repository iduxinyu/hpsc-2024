#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }


  for(int i=0; i<N; i++) {
  /* for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
   */

    __m256 xivec=_mm256_set1_ps(x[i]);
     __m256 yivec=_mm256_set1_ps(y[i]);
    
     __m256 xjvec=_mm256_load_ps(x);
     __m256 yjvec=_mm256_load_ps(y);
    
     __m256 mvec=_mm256_load_ps(m);
    
    // __m256 ivec=_mm256_set1_ps(i);
    // __m256 jvec=_mm256_set_ps(0,1,2,3,4,5,6,7);
    // __mmask8 mask=_mm256_cmp_ps_mask(ivec,jvec,0);
    
     __m256 zero=_mm256_set1_ps(0);

     __m256 rxvec=_mm256_sub_ps(xivec,xjvec);
     __m256 ryvec=_mm256_sub_ps(yivec,yjvec);
    
     __m256 rx2vec=_mm256_mul_ps(rxvec,rxvec);
     __m256 ry2vec=_mm256_mul_ps(ryvec,ryvec);

     __m256 rvec=_mm256_add_ps(rx2vec,ry2vec);

     __mmask8 mask=_mm256_cmp_ps_mask(rvec,zero,_MM_CMPINT_NE);

     rvec=_mm256_rsqrt14_ps(rvec);
     
     __m256 r2vec=_mm256_mul_ps(rvec,rvec);
     __m256 r3vec=_mm256_mul_ps(r2vec,rvec);

     

     __m256 fxsub=_mm256_mul_ps(rxvec,mvec);
     fxsub=_mm256_mul_ps(fxsub,r3vec);
     fxsub=_mm256_mask_blend_ps(mask, zero, fxsub);
     __m512 fxsub512=_mm512_castps256_ps512(fxsub);
     fx[i] -= _mm512_reduce_add_ps(fxsub512);


     __m256 fysub=_mm256_mul_ps(ryvec,mvec);
     fysub=_mm256_mul_ps(fysub,r3vec);
     fysub=_mm256_mask_blend_ps(mask, zero, fysub);
     __m512 fysub512=_mm512_castps256_ps512(fysub);
      fy[i] -= _mm512_reduce_add_ps(fysub512);
   


 //    printf("fxisub %g  fyisub %g \n",fxisub,fyisub);
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
