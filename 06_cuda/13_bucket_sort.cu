#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucketSort(int* key, int* bucket, int range){

    int id=blockIdx.x*blockDim.x+threadIdx.x;

    

    if(id<range){
	bucket[id]=0;
    }

    __syncthreads();
    atomicAdd(&bucket[key[id]],1);

    int offset0=0, offset1=0;
    for(int i=0;i<range;i++){

	offset0=offset1;    
	offset1+=bucket[i];
	
	if(id<offset1 && id>=offset0){
	key[id]=i;

	
	continue;
	
	}
    }

}

int main() {
  int n = 50;
  int range = 5;

  int* key;
  int* bucket;
  cudaMallocManaged(&key,n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));

//  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
/*
  std::vector<int> bucket(range); 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
*/
  bucketSort<<<1,n>>>(key,bucket,range);
  cudaDeviceSynchronize();
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
