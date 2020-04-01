#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define RAND01 ((double)random() / (double)RAND_MAX)

//Random Function
void random_fill_LR(double **matrix1, double **matrix2, int nU, int nI, int nF)
{
  srandom(0);

  for(int i = 0; i < nU; i++)
    for(int j = 0; j < nF; j++)
      matrix1[i][j] = RAND01 / (double) nF;

  for(int i = 0; i < nF; i++)
    for(int j = 0; j < nI; j++)
      matrix2[i][j] = RAND01 / (double) nF; 
  
}

// N1 -> rows of matrix1
// N2 -> columns of matrix2
// M -> columns of matrix1 and rows of matrix2
void matrix_product(double **matrix1, double **matrix2, double **result_matrix, int N1, int N2, int M){
    int k; // rows and columns of result_matrix
    //determina la matrice prodotto
    for(int i=0; i<N1; i++){ //i: blocca la riga di m1 da moltiplicare m1[i][?]
        for(int j=0; j<N2; j++){ //j: blocca la colonna di m2 da moltiplicare m2[?][j]
            result_matrix[i][j] = 0;
                for(k=0; k<M; k++){
                    result_matrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }
        }
    }
}

void copy_matrix(double **matrix,double **original_matrix,int rows,int cols){
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
        matrix[i][j] = original_matrix[i][j];
}

//rows->rows of non-zero elements of A
//cols->cols of non-zero elements of A
//nz-> number of non-zero elements
void matrix_factorization(double **L,double **R, double **B,double **A,double **matrix1, double **matrix2,int *rows_nz,int *cols_nz,int nz,float alpha,int nF,int nU,int nI){

    int j=0;
    for(int i=0; i<nz; i++){
        for(int k=0; k<nF;k++){
            L[rows_nz[i]][k]= L[rows_nz[i]][k]-2*alpha*(A[rows_nz[i]][cols_nz[j]]-B[rows_nz[i]][cols_nz[j]])*(-matrix2[k][cols_nz[j]]);
            R[k][cols_nz[j]]= R[k][cols_nz[j]]-2*alpha*(A[rows_nz[i]][cols_nz[j]]-B[rows_nz[i]][cols_nz[j]])*(-matrix1[rows_nz[i]][k]);
        }
        j++;
    }
 
  copy_matrix(matrix1,L,nU,nF);
  copy_matrix(matrix2,R,nF,nI);
}

//Main
int main(int argc, char *argv[]){

  // check file
  if (argv[1]==NULL){
    std::cout<<"Error: insert file!" << std::endl;
    return -1;
  }

  std::ifstream OpenFile(argv[1]);

  std::string line;
  double data[3]; //vector for first three parameter

  // Initialization matrix L, R, A
  double **L;
  double **R;
  double **A;
  double **B;

  //Matrix initialization
  double **L_original;
  double **R_original;

  int iter; // number of iterations
  float alpha; // value of convergence rate
  int features;
  int rows,cols; // values of matrix A
  int no_zero; // number of non-zero elements into matrix A
  int rows_A, cols_A; // values of max matrix A
  double value_A; // value of elements into matrix A
  int *rows_nz; // value row no-zero
  int *cols_nz; // value column no-zero
  double diff_pow; // value of pow differences no-zero
  double sum_diff=0; // value sum of pow_diff

  int while_i=0; // value line of file
  int while_j=0; // increment value of no-zero element

  while (std::getline(OpenFile, line))
  {
      std::istringstream iss(line);
      if(while_i<3)

        iss >> data[while_i]; // save value of the first 3 line of file

      else if(while_i==3){

        if (!(iss >> rows_A >> cols_A >> no_zero)) { break; } //error

            // Inizialitation of matrix A
            A = new double*[rows_A]();
            for(int i=0; i<rows_A; i++)
                A[i] = new double[cols_A]();

            // inizialitation value rows and columns no-zero
            rows_nz = new int[no_zero];
            cols_nz = new int[no_zero];

      }else{

        iss >> rows >> cols >> value_A;            

        // Insert values into matrix A
        A[rows][cols]={value_A};

        // insert values rows and columns no-zero
        rows_nz[while_j] = rows;
        cols_nz[while_j] = cols;

        while_j++;

      }

      while_i++;
  }

  // Insert value into variables iteration, alpha, features
  iter = data[0];
  alpha = data[1];
  features = data[2];

  // Values of rows and columns of matrix L_original and R_original
  int nU = rows_A;
  int nI = cols_A;
  int nF = features;

  // Matrix original L initialization
  L_original = new double*[nU]();
    for(int i=0; i<nU; i++)
        L_original[i] = new double[nF]();

  // Matrix original R initialization
  R_original = new double*[nF]();
    for(int i=0; i<nF; i++)
        R_original[i] = new double[nI]();

  // Matrix L initialization
  L = new double*[nU]();
    for(int i=0; i<nU; i++)
        L[i] = new double[nF]();

  // Matrix R initialization
  R = new double*[nF]();
    for(int i=0; i<nF; i++)
        R[i] = new double[nI]();

    //Matrix B initialization
  B = new double*[nU]();
    for(int i=0; i<nU; i++)
        B[i] = new double[nI]();

  // Fill matrix L_origianl and R_original with random value
  random_fill_LR(L_original, R_original, nU, nI, nF);

  // Copy matrix L_original into L
  copy_matrix(L,L_original,nU,nF);
  // Copy matrix R_original into R
  copy_matrix(R,R_original,nF,nI);

  matrix_product(L, R, B, nU, nI, nF);
  
  for(int i=0;i<iter;i++){

    matrix_factorization(L,R,B,A,L_original,R_original,rows_nz,cols_nz,no_zero,alpha,nF,nU,nI);
    matrix_product(L, R, B, nU, nI, nF);
  }

  double max = 0;
  double foo = 0;
  int col = 0;

  for(int i=0; i<nU; i++){
      for(int j=0; j<nI; j++){
          if(A[i][j] == 0){

            // Search max into B matrix
            foo = B[i][j];

            if(max < foo){
                max = foo;
                col = j;
            }
          }
      }
      max = 0;

      std::cout<<col<<std::endl;
     
  }

  return 0;
}
