
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <chrono>
#include <math.h>

using namespace std;
typedef vector<vector<float>> matrix;


__global__ void init(float *u, float *v,float *p,float *b,float *un, float *vn, float *pn)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;

	//int j=id/nx;
	//int i=id%nx;

	u[id]=0;
	v[id]=0;
	p[id]=0;
	b[id]=0;
	un[id]=0;
	vn[id]=0;
	pn[id]=0;
}

__global__ void comB(float *u, float *v,float *b, int nx, int ny, double rho, double dx, double dy, double dt)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=1&&j<ny-1&&i>=1&&i<nx-1)
	{
		b[j*nx+i]=rho * (1 / dt *
					((u[j*nx+i + 1] - u[j*nx+i - 1]) / (2 * dx) + (v[(j + 1)*nx+i] - v[(j - 1)*nx+i]) / (2 * dy)) -
					pow((u[j*nx+i + 1] - u[j*nx+i - 1]) / (2 * dx), 2) - 2 * ((u[(j + 1)*nx+i] - u[(j - 1)*nx+i]) / (2 * dy) *
					(v[j*nx+i + 1] - v[j*nx+i - 1]) / (2 * dx)) - pow((v[(j + 1)*nx+i] - v[(j - 1)*nx+i]) / (2 * dy), 2));
	}

}

__global__ void comP(float *p, float *pn, float *b, int nx, int ny, double dx, double dy, double dt)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=1&&j<ny-1&&i>=1&&i<nx-1)
	{
		p[id] = (dy * dy * (pn[j*nx+i + 1] + pn[j*nx+i - 1]) +
						dx * dx * (pn[(j + 1)*nx+i] + pn[(j - 1)*nx+i]) -
						b[j*nx+i] * pow(dx, 2) * pow(dy, 2))
						/ (2 * (pow(dx, 2) + pow(dy, 2)));
	}

}

//可以并行吗？？ 有依赖性吗???????
__global__ void boundryP(float *p, int ny, int nx)
{

	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=0&&j<ny)
	{
		p[j*nx+nx - 1] = p[j*nx+nx - 2];
		p[j*nx] = p[j*nx+1];
	}

	if(i>=0&&i<nx)
	{
		p[i] = p[nx+i];
		p[(ny - 1)*nx+i] = 0;
	}
}

__global__ void copyP(float* p, float* pn, int ny, int nx)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	pn[id]=p[id];

}

__global__ void copyUV(float *u, float *v, float *un, float *vn, int ny, int nx)
{

	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=0&&j<ny&&i>=0&&i<nx)
	{
		un[id]=u[id];
		vn[id]=v[id];
	}
	
}

__global__ void comUV(float *u,float *v, float *un,float *vn, float *p, double dx,double dy, double dt, double rho, double nu,int nx, int ny)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=1&&j<ny-1&&i>=1&&i<nx-1)
	{
		u[id] = un[id] - un[id] * dt / dx * (un[id] - un[id - 1])
					- un[id] * dt / dy * (un[id] - un[(j - 1)*nx+i])
					- dt / (2 * rho * dx) * (p[(j + 1)*nx+i] - p[(j - 1)*nx+i])
					+ nu * dt / pow(dx, 2) * (un[id + 1] - 2 * un[id] + un[id - 1])
					+ nu * dt / pow(dy, 2) * (un[(j + 1)+i] - 2 * un[id] + un[(j - 1)*nx+i]);

		
		v[id] = vn[id] - vn[id] * dt / dx * (vn[id] - vn[id - 1])
					- vn[id] * dt / dy * (vn[id] - vn[(j - 1)*nx+i])
					- dt / (2 * rho * dx) * (p[(j + 1)*nx+i] - p[(j - 1)*nx+i])
					+ nu * dt / pow(dx, 2) * (vn[id + 1] - 2 * vn[id] + vn[id - 1])
					+ nu * dt / pow(dy, 2) * (vn[(j + 1)*nx+i] - 2 * vn[id] + vn[(j - 1)*nx+i]);
	}

}

__global__ void boundryUV(float *u, float *v, int nx, int ny)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;

	int j=id/nx;
	int i=id%nx;

	if(j>=0&&j<ny)
	{
		u[j*nx] = 0;
		u[j*nx+nx - 1] = 0;
		v[j*nx] = 0;
		v[j*nx+nx - 1] = 0;
	}

	
	if(i>=0&&i<nx)
	{
			u[i] = 0;
			u[(ny - 1)*nx+i] = 1;
			v[i] = 0;
			v[(ny - 1)*nx+i] = 0;
	}

	u[id]=id;
	v[id]=id;

}

int main() {
	int nx = 41;
	int ny = 41;
	int nt = 500;
	int nit = 50;
	double dx = 2. / (nx - 1);
	double dy = 2. / (ny - 1);
	double dt = .01;
	double rho = 1.;
	double nu = .02;

	float *u,*v,*p,*b,*un,*vn,*pn;

	cudaMallocManaged(&u,ny*nx*sizeof(float));
	cudaMallocManaged(&v,ny*nx*sizeof(float));
	cudaMallocManaged(&p,ny*nx*sizeof(float));
	cudaMallocManaged(&b,ny*nx*sizeof(float));
	cudaMallocManaged(&un,ny*nx*sizeof(float));
	cudaMallocManaged(&vn,ny*nx*sizeof(float));
	cudaMallocManaged(&pn,ny*nx*sizeof(float));

	init<<<ny,nx>>>(u, v, p , b , un, vn, pn);
	cudaDeviceSynchronize();

	ofstream ufile("u.dat");
	ofstream vfile("v.dat");
	ofstream pfile("p.dat");


	for (int n = 0; n < nt; n++) {
		
		comB<<<ny,nx>>>(u, v, b, nx, ny, rho, dx, dy, dt);
		cudaDeviceSynchronize();


		for (int it = 0; it < nit; it++) {
			copyP<<<ny,nx>>>(p, pn, ny, nx);
			cudaDeviceSynchronize();

			comP<<<ny,nx>>>(p, pn, b , nx, ny, dx, dy, dt);
			cudaDeviceSynchronize();

			boundryP<<<ny,nx>>>(p, ny, nx);
			cudaDeviceSynchronize();
			
		}
		
		copyUV<<<ny,nx>>>(u, v, un, vn, ny, nx);
		cudaDeviceSynchronize();



		comUV<<<ny,nx>>>(u, v, un, vn,p, dx, dy, dt,  rho, nu, nx, ny);
		cudaDeviceSynchronize();

		
		printf("n =%d \n", n);

		boundryUV<<<ny,nx>>>(u, v, nx, ny);
		cudaDeviceSynchronize();

		if (n % 10 == 0) {
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
					ufile << u[j*nx+i] << " ";
			ufile << "\n";
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
					vfile << v[j*nx+i] << " ";
			vfile << "\n";
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
					pfile << p[j*nx+i] << " ";
			pfile << "\n";
		}
	}



	ufile.close();
	vfile.close();
	pfile.close();
}
