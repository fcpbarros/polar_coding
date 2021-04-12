/* Polar Codes
 * Sucessive Cancellation Decoding: transmission over AWGN or Rayleigh
 * Date: Apr/15/2019
 * Author: Waslon T. A. Lopes           */

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "rand.h"
using namespace std;

// Function prototypes

int **alloc_nodes(int N);
void free_nodes(int **nodes, int N);
double **aloca_beliefs(int N);
void free_beliefs(double **beliefs, int N);
int **aloca_ucap(int N);
void free_ucap(int **ucap, int N);
double f_function(double r1, double r2);
double g_function(double r1, double r2, int b);
void intercambio(double& x, double& y);
void intercambio2(int& x, int& y);
void ordenar (double [], int [], int);
int *channel_reliability(int N, double e);

// Main program
int main(){

int *A, *B,*F, Rayleigh;
double alpha, *R, No, EbNodB, EbNo, EbNo_min, EbNo_max, EbNo_step, rate, sigma, e;
uint64_t  nsim_min, nsim_max, MaxErrors;
int N, K;
char filename[60];

reset_rand(); //reset the random number generator

cout << "Code dimension (N)? ";
cin  >> N;
cout << "Number of information bits (K)? ";
cin  >> K;
cout << "Minimum Eb/No (EbNo_min)? ";
cin  >> EbNo_min;
cout << "Maximum Eb/No (EbNo_max)? ";
cin  >> EbNo_max;
cout << "Increment Eb/No (EbNo_step)? ";
cin  >> EbNo_step;
cout << "AWGN (0) or Rayleigh (1)? ";
cin  >> Rayleigh;
cout << "Number of errors for calculation of BER (MaxErrors)? ";
cin  >> MaxErrors;
cout << "Number (minimum) of simulations (nsim_min)? ";
cin  >> nsim_min;
cout << "Number (maximum) of simulations (nsim_max)? ";
cin  >> nsim_max;
cout << "Filename for data recording? ";
cin  >> filename;
cout << endl;
 


// Test if N is power of 2
if((N&(N-1))!=0) {cout<< "N should be power of two!\n"; exit(1);}

// File opening
std::ofstream output_file(filename,ios::app);
if(!output_file)
{std::cout << "The output file can not be opened!";
return 1;}

int *data = new int[K]; //Length of data is K
int *received_data = new int[K]; //Length of data is K

A = new (nothrow) int[N]; if (!A){cout << "Memory allocation failed\n";}
B = new (nothrow) int[N]; if (!B){cout << "Memory allocation failed\n";}
F = new (nothrow) int[N]; if (!F){cout << "Memory allocation failed\n";} //Frozen bits positions
R = new (nothrow) double[N]; if (!R){cout << "Memory allocation failed\n";}
//alpha = new (nothrow) double[N]; if (!alpha){cout << "Memory allocation failed\n";} //Fading gain

double **beliefs = aloca_beliefs(N); 
int **ucap = aloca_ucap(N); 
int **nodes = alloc_nodes(N);
int depth = (int)((float)log(N)/log(2))+1; //Tree depth

e = exp(-((double)K/N)* pow(10,2.0/10)); // Battacharaya
//e = ((double)K/N)* pow(10,2.0/10);

//Frozen bits positions 
int* order = channel_reliability(N, e);

for(int i=0; i<(N-K); i++)
  F[order[i]] = 1; // F[i] = 1 (frozen)


for(EbNodB = EbNo_min; EbNodB <= EbNo_max ; EbNodB += EbNo_step){
//EbNodB = 3.0;
rate = (double)K/N;
EbNo =  pow(10,EbNodB/10) ;
sigma = sqrt(1.0/(EbNo*N/K)); //Rate here 


uint64_t  NumberOfBlocks = 0;
double    SumOfMedia = 0.0;
uint64_t  TotalErrors = 0;
uint64_t  errors = 0;


while(!(((NumberOfBlocks >= nsim_min)&&(TotalErrors >= MaxErrors))||(NumberOfBlocks>=nsim_max))){
  //Load the random bits
  
  for(int i=0;i<K;i++){
		data[i] = gera_bits(0.5);
  }

  //Load vector A with random bits + Frozen bits
  int pos=0;
  for(int i=0;i<N;i++){
		if(F[i])
		  A[i] = 0;
		else{
		  A[i] = data[pos];
		  pos++;
		}
	 }


  //Polar encoding of vector A

  int steps = (int) ((float)log(N)/log(2));
  for(int L=0;L<steps;L++){
	 int LE = (int)pow(2,L+1);
	 int LE1 = LE/2;
	 for(int J=1;J<(LE1+1);J++){
		 for(int I=J;I<(N+1);I+=LE){
			  int IP = I + LE1;
			  A[I-1] = (A[IP-1] + A[I-1])%2;
		  }
	 }
  } 

  //Transmission of Vector A over Channel (AWGN or Rayleigh)

  if(Rayleigh){
	 for(int i=0;i<N;i++){
	 alpha = sqrt(pow(gauss(0,0.5),2)+pow(gauss(0,0.5),2));
		  R[i] = alpha*(1-2*A[i]) + gauss(0,pow(sigma,2));
	 }
  }
  else{
	 for(int i=0;i<N;i++)
		  R[i] = (1-2*A[i]) + gauss(0,pow(sigma,2));
  }

  //Decoding
   
  if(Rayleigh){
	 for(int i=0; i<N; i++)
		//beliefs[0][i]=R[i]; // Copy the channel outputs
		beliefs[0][i]=2*R[i]/(alpha*pow(sigma,2.0)); // Copy the channel outputs
  }
  else{
	 for(int i=0; i<N; i++)
		beliefs[0][i]=2*R[i]/pow(sigma,2.0); // Copy the channel outputs
  }

  int temp;
  int shift;

  int done = 0; // Initialization
  int posX, posY; // Position of the node
  posX = 0; posY = 0; // Start at root node
  for(int i = 0; i < depth; i++)
	 for(int j =0; j < pow(2,i); j++)
	 //for(int j =0; j < N; j++)
		nodes[i][j]=0;

  while(done == 0){
	 if(posX == depth-1){ //Node or LEAF? It is a Leaf
		// Make decision
		if(F[posY]) //Frozen
		  ucap[posX][posY] = 0;
		else{
		  if (beliefs[posX][posY]>=0)
			 ucap[posX][posY] = 0;
		  else
			 ucap[posX][posY] = 1;
		  }

		nodes[posX][posY]=3; //State update

		if(posY==(int)pow(2,depth-1)-1){
		  done = 1;
		}
		else{
		  posX--; posY /= 2; //Next node
		  }
		}
	 else //{//NODE or leaf? It is a node
		switch(nodes[posX][posY]){
		  case 0: //Goes to left leaf
			 //Processing of beliefs
			 temp = N/pow(2,posX+1);
			 shift =  N/pow(2,posX) * posY;
			 for(int i=0; i<temp ; i++)
				  beliefs[posX+1][i+shift] = f_function(beliefs[posX][i+shift],beliefs[posX][i+temp+shift]);
			 nodes[posX][posY]=1; //State update
			 posX++; posY *= 2;   //Next node
			 break;
		  case 1: //Goes to right leaf
			 //Processing of beliefs
			 temp = N/pow(2,posX+1);
			 shift =  N/pow(2,posX) * posY;
			 for(int i=0; i<temp ; i++)
				beliefs[posX+1][i+shift+temp] = g_function(beliefs[posX][i+shift],beliefs[posX][i+temp+shift],ucap[posX+1][i+shift]);
			 nodes[posX][posY]=2;     //State update
			 posX++; posY = posY*2+1; //Next node
			 break;
		  case 2: //Goes up
			 //Sent bits up
			 temp = N/pow(2,posX);
			 shift =  temp * posY;
			 for(int i=0; i<temp/2; i++)
				ucap[posX][i+shift] =(ucap[posX+1][i+shift]+ucap[posX+1][i+temp/2+shift])%2; 
			 for(int i=temp/2; i<temp; i++)
				ucap[posX][i+shift] = ucap[posX+1][i+shift]; 
			 nodes[posX][posY]=3; //State update
			 posX--; posY /= 2;   //Next node
			 break;
		} 
  } // End of decoding
  
  //Extract  Frozen bits from the decoded vector

    int pos1=0, pos2=0;
		for(int i=0;i<N;i++){
		  if(F[i])
			 pos1++;
		  else{
			 received_data[pos2] = ucap[depth-1][pos1];
			 pos1++;
			 pos2++;
		  }
		}


  // Errors update

  errors=0;
  for(int i = 0; i < K; i++)
	 if(received_data[i]!=data[i]){
		  errors++;
		  TotalErrors++;
	 }


// for(int i = 0; i < N; i++)
//    if(!F[i])
//		if(B[i]!=ucap[depth-1][i]){
//		  errors++;
//		  TotalErrors++;
//		}

	NumberOfBlocks++;

	SumOfMedia += (1.0*errors)/(1.0*K);


} //End of nsim


//cout << EbNodB << "\t" << errors/(double(nsim*N*rate)) << endl;
//output_file << EbNodB << "\t" << errors/(double(nsim*N*rate)) << endl;

cout << EbNodB << "\t" << SumOfMedia/NumberOfBlocks << endl;
output_file << EbNodB << "\t" << SumOfMedia/NumberOfBlocks << endl;

} //End of EbNodB

//Free memory
delete[] A;
delete[] B;
delete[] R;
delete[] F;
delete[] order;
free_nodes(nodes,N);
free_beliefs(beliefs,N);
free_ucap(ucap,N);

return 0;
} // End of main()


int **alloc_nodes(int N){
  int depth = (int)(log(N)/log(2))+1;
  int **nodes = new int*[depth];
  for (int i=0; i < depth; i++)
	  nodes[i] = new int[(int)pow(2,i)];
  return nodes;
};


void free_nodes(int **nodes, int N){
  int depth = (int)(log(N)/log(2))+1;
  for (int i=0; i < depth; i++)
	  delete [] nodes[i];
  delete [] nodes;
  return ;
};

double **aloca_beliefs(int N){
  int depth = (int)(log(N)/log(2))+1;
  double **beliefs = new double*[depth];
  for (int i=0; i < depth; i++)
	  beliefs[i] = new double[N];
  return beliefs;
};

void free_beliefs(double **beliefs, int N){
  int depth = (int)(log(N)/log(2))+1;
  for (int i=0; i < depth; i++)
	  delete [] beliefs[i];
  delete [] beliefs;
  return;
};

int **aloca_ucap(int N){
  int depth = (int)(log(N)/log(2))+1;
  int **ucap = new int*[depth];
  for (int i=0; i < depth; i++)
	  ucap[i] = new int[N];
  return ucap;
};

void free_ucap(int **ucap, int N){
  int depth = (int)(log(N)/log(2))+1;
  for (int i=0; i < depth; i++)
	  delete [] ucap[i];
  delete [] ucap;
  return;
};

double f_function(double r1, double r2){
  double signal_r1, signal_r2, minor;

  if(r1>0)
	 signal_r1 = 1;
  else
	 signal_r1 = -1;

  if(r2>0)
	 signal_r2 = 1;
  else
	 signal_r2 = -1;

  if(abs(r1) < abs(r2))
	 minor = abs(r1);
  else
	 minor = abs(r2);
	 
  return (signal_r1 * signal_r2 * minor);
}

double g_function(double r1, double r2, int b){
  return (r2 + (1-2.0*b)*r1);
}

void intercambio(double& x, double& y){
double temp = x;
x = y;
y = temp;
}

void intercambio2(int& x, int& y){
int temp = x;
x = y;
y = temp;
}

void ordenar(double a[], int c[], int n){
for(int i = n; i > 0; i--)
  for(int j = 0; j < i-1 ; j++)
	 if(a[j] > a[j+1]){
		intercambio(a[j], a[j+1]);
		intercambio2(c[j], c[j+1]);
	 }
}

int *channel_reliability(int N, double e){
// This function returns a int vector with the indices of subchannels ordered in
// accordance with the reliability of each channell

int steps = (int) (log(N)/log(2)); // 2^n = N

double *capacity = new (nothrow) double[N];
if (!capacity){cout << "Memory allocation failed\n";}

int *order = new (nothrow) int[N];
if (!order){cout << "Memory allocation failed\n";}

for(int i=0;i<N;i++)
  order[i]=N-i-1;

capacity[0] = 2*e - pow(e,2);
capacity[1] = pow(e,2);

for(int i=1;i<steps;i++)
  for(int j=pow(2,i);j>0;j--){
	 capacity[2*(j-1)+1] = pow(capacity[j-1],2);
	 capacity[2*(j-1)] = 2*capacity[j-1] - pow(capacity[j-1],2);
  }


ordenar(capacity, order, N);

delete[] capacity;

return(order);
}

