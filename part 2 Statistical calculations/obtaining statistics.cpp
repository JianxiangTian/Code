//CPP file for obtaining statistics of the 2D configuration, including:
//     (1) pair correlation functions
//     (2) number variance
//     (3) structure factor 

using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <random>


#define N 10000 //Number of particles
#define SD 2 //dimension
#define Bin 0.2 //shell thickness

double Center[N][SD];
double Lx;  //side length of box
double Ly;
int Nbin; //for g(r)

int Nk = 600;
double Kbin = 0.08;
#define Kbin_num 200

double pi = 3.14159265358979;
#define Rbin 0.2  //for number variance
#define Nr 100000 //sampling points for number variance at each r
#define MAXY 32767 

void read_config()
{
	FILE * fp;
	if ((fp = fopen("./vertex_relax.txt", "r")) == NULL)
	{
		perror("Cannot open file!\n");
		exit(1);
	}
	else
	{
		double temp_t;
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Lx);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &temp_t);
		fscanf(fp, "%lf", &Ly);
		cout << "Lx: " << Lx << endl;
		cout << "Ly: " << Ly << endl;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < SD; j++)
			{
				fscanf(fp, "%lf", &Center[i][j]);
			}

		Nbin = (int)floor(Ly / 2.0 / Bin);
		fclose(fp);
	}
}



double MinDis(int m, int n)
{
  //find the minimal distance between the centers of two polyhedra in Ecludean space...
  //by checking all the images of Poly[m], while keeping Poly[n] in the central box
  //record the index of the box which Poly[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double dx = Center[m][0] - Center[n][0];
  double dy = Center[m][1] - Center[n][1];
 

  double dist = 1000000000.0; //just a large number...


  //loop over all possible images of Point m, keep n fixed in the center simulation box....
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
	{
	  double tempd[SD]; 
	  tempd[0] = dx + (double)i*Lx;
	  tempd[1] = dy + (double)j*Ly;
	  
	  double tempdist = tempd[0]*tempd[0]+tempd[1]*tempd[1]; 
	  
	  
	  if(tempdist < dist) // store the smallest distance...
	    {
			dist = tempdist; 
	    }
	  
	}
  
  
  return sqrt(dist); //this is center-to-center distance...
}

double Get_NDensity() // the number density, for computing g2...
{
  double VLambda;
  VLambda = Lx * Ly;

  return (double)N/VLambda;
}

void Get_PairCorr()
{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //compute the minimum length of the lattice vectors
  double * g;
  
  g = new double [Nbin];
  for(int i = 0; i < Nbin; i++)
  {
		g[i] = 0.0;
   }	
  printf("Computing G2 and g2 now....\n");

  
  //loop over all particles
  for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
	{
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the radial g2....
	  double temp_dis = MinDis(j, i)/Bin;
	  int int_dis = (int)floor(temp_dis);
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(j != i && temp_dis <= Nbin)
	    {
	      
	      //double delta_dis = temp_dis - int_dis;
	      //if(delta_dis>0.5) int_dis++;
	      
	      g[int_dis] = g[int_dis] + 1.0;
	    }
	  }

	 
  double rho_n = Get_NDensity();

  for(int r = 0; r < Nbin; r++)
    {
      g[r] = g[r] / (N * rho_n * pi * ((r+1.0) * (r+1.0) - r * r) * Bin * Bin);
    }

 
  
  FILE* fp = fopen("./g2_relax.txt","w");
  for(int r = 0; r < Nbin; r++)
    fprintf(fp, "%lf\t%lf\n", (r + 0.5) * Bin, g[r]);
  fclose(fp);
  delete [] g;
}


void num_var()
{
	printf("Computing number variance now....\n");
	int Ns = (int)floor(Ly / 4.0 / Rbin) - 1;
	cout << "Ns = " << Ns << endl;
	double* sigma = new double [Ns];
	double ave, square_ave;
	int Nc;
	double r = Rbin;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(0.0, 1.0);
	
	for(int t = 0; t < Ns; t++)
	{
		cout << "r = " << r << endl;
		ave = 0.0;
		square_ave = 0.0;
		for(int n = 0; n < Nr; n++)
		{
			Nc = 0;
			double cx = distribution(gen) * Lx;
			double cy = distribution(gen) * Ly;
			for(int i = 0; i < N; i++)
			{
				double dx = fabs(Center[i][0] - cx);
				if (dx > Lx / 2)
				{
					dx = Lx - dx;
				}
				double dy = fabs(Center[i][1] - cy);
				if (dy > Ly / 2)
				{
					dy = Ly - dy;
				}
				double dist = sqrt(dx * dx + dy * dy);
				if(dist < r)
				{
					Nc++; 
				}
			}
			ave = ave + (double)Nc;
			square_ave = square_ave + (double)(Nc * Nc);
		}
		ave = ave / Nr;
		square_ave = square_ave / Nr;
		cout << "ave = " << ave << endl;
		cout << "square_ave = " << square_ave << endl;
		sigma[t] = square_ave - ave * ave;
		r = r + Rbin;
	}
	FILE* fp = fopen("./nv_relax.txt","w");
	for(int r = 0; r < Ns; r++)
		fprintf(fp, "%lf\t%lf\n", (r + 1) * Rbin, sigma[r]);
	fclose(fp);
	delete [] sigma;
}

double GetInnerProduct(double Vector1[SD], double Vector2[SD])
{
	double sum = 0;

	for(int i = 0; i < SD; i++)
		sum += Vector1[i] * Vector2[i];
	return sum;
}

double Get_Sk(double Vector_k[SD])
{
	double Sk;
	double sum_cos = 0.0;
	double sum_sin = 0.0;
	double temp;
	for(int i = 0; i < N; i++)
	{
			temp = GetInnerProduct(Center[i], Vector_k); 
			sum_cos += cos(temp);
			sum_sin += sin(temp);
	}
	Sk = (sum_cos * sum_cos + sum_sin * sum_sin) / (double)N;
	return Sk;
}

void Print_Sk(double K_Histo[Kbin_num], int K_Counter[Kbin_num])
{
	ofstream histo_out;
	histo_out.open("./Sk_relax.txt");
	for(int t = 0; t < Kbin_num; t++)
		if(K_Counter[t] > 0)
			histo_out << (t + 0.5) * Kbin << " " << K_Histo[t] << endl;
	histo_out.close();
}

void Get_KHistogram()
{
	printf("Computing Sk now....\n");
	double Sk, k_dis;
	int t;
	double KPoint[SD];
	double K_Histo[Kbin_num];
	int K_Counter[Kbin_num];
	for(int i = 0; i < Kbin_num; i++)
	{
		K_Histo[i] = 0.0;
		K_Counter[i] = 0;
	}
	for(int i = 1; i <= Nk; i++)
		for(int j = - Nk; j <= Nk; j++)
		{
			KPoint[0] = i * 2.0 * pi / Lx;
			KPoint[1] = j * 2.0 * pi / Ly;
			Sk = Get_Sk(KPoint);
			k_dis = 0.0;
			for(int d = 0; d < SD; d++)
				k_dis += KPoint[d] * KPoint[d];
			k_dis = sqrt(k_dis);
			t = floor(k_dis / Kbin);
			if(t < Kbin_num)
			{
				K_Histo[t] += Sk;
				K_Counter[t] ++;
			}
		}
	for(int i = 1; i <= Nk; i++)
	{
		KPoint[0] = 0.0;
		KPoint[1] = i * 2.0 * pi / Ly;
		Sk = Get_Sk(KPoint);
		k_dis = 0.0;
		k_dis = KPoint[1];
		t = floor(k_dis / Kbin);
		if(t < Kbin_num)
		{
			K_Histo[t] += Sk;
			K_Counter[t] ++;
		}
	}
	for(t = 0; t < Kbin_num; t++)
		if(K_Counter[t] != 0)		
			K_Histo[t] = K_Histo[t] / K_Counter[t];
	Print_Sk(K_Histo, K_Counter);
}

int main()
{
	read_config();	
	Get_PairCorr();
	Get_KHistogram();
	srand(time(NULL));
	num_var();
	return 1;
}

