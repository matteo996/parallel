#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define RAND01 ((double)random() / (double)RAND_MAX)

//Matrix initialization
double **L_original;
double **R_original;

//Random Function
void random_fill_LR(int nU, int nI, int nF)
{
  srandom(0);
  for(int i = 0; i < nU; i++)
    for(int j = 0; j < nF; j++)
      L_original[i][j] = RAND01 / (double) nF;

  for(int i = 0; i < nF; i++)
    for(int j = 0; j < nI; j++)
      R_original[i][j] = RAND01 / (double) nF; // CORREGGERE CON COPIA
}

//nr-> number of rows
//nc-> number of columns
void print_matrix(double **matrix, int nr, int nc){
  std::cout<<" ("<<nr<<"x"<<nc<<")"<<std::endl;
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++)
      std::cout<<matrix[i][j]<<"\t";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
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
void matrix_factorization(double **L,double **R, double **B,double **A,int *rows_nz,int *cols_nz,int nz,float alpha,int nF,int nU,int nI){
  int j=0;
  for(int i=0; i<nz; i++){
      for(int k=0; k<nF;k++){
        L[rows_nz[i]][k]= L[rows_nz[i]][k]-2*alpha*(A[rows_nz[i]][cols_nz[j]]-B[rows_nz[i]][cols_nz[j]])*(-R_original[k][cols_nz[j]]);
        R[k][cols_nz[j]]= R[k][cols_nz[j]]-2*alpha*(A[rows_nz[i]][cols_nz[j]]-B[rows_nz[i]][cols_nz[j]])*(-L_original[rows_nz[i]][k]);
      }
    j++;
  }
  copy_matrix(L_original,L,nU,nF);
  copy_matrix(R_original,R,nF,nI);
}

//Main
int main(int argc, char *argv[]){

  // check file
  if (argv[1]==NULL){
    std::cout<<"Error: insert file!" << std::endl;
    return -1;
  }

  else if (argc > 2){
    std::cout<< "Error: wrong number of arguments!" << std::endl;
    return -1;
  }

  //std::ifstream OpenFile("InputData.txt");
  std::ifstream OpenFile(argv[1]);

  // take the name of file and remove extension
  std::string name_file = argv[1];
  size_t lastindex = name_file.find_last_of(".");
  name_file = name_file.substr(0, lastindex);

  std::string line;
  double data[3]; //vector for first three parameter

  // Initialization matrix L, R, A
  double **L;
  double **R;
  double **A;

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

  std::ofstream myfile; // variable file where write final value

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

  // Print matrix A
  // std::cout<<"Matrix A:";
  // print_matrix(A,rows_A,cols_A);

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

  // Fill matrix L_origianl and R_original with random value
  random_fill_LR(nU,nI,nF);

  // Copy matrix L_original into L
  copy_matrix(L,L_original,nU,nF);
  // Copy matrix R_original into R
  copy_matrix(R,R_original,nF,nI);

  // Print matrix L_original and R_original
  // std::cout<<"Matrix L original:";
  // print_matrix(L_original,nU,nF);
  // std::cout<<"Matrix R original:";
  // print_matrix(R_original,nF,nI);

  //Matrix B initialization
  double **B = new double*[nU]();
    for(int i=0; i<nU; i++)
        B[i] = new double[nI]();

  matrix_product(L, R, B, nU, nI, nF);

  // std::cout<<"Matrix B original:";
  // print_matrix(B,nU,nI);

  for(int i=0;i<iter;i++){

    matrix_factorization(L,R,B,A,rows_nz,cols_nz,no_zero,alpha,nF,nU,nI);
    matrix_product(L, R, B, nU, nI, nF);

  }

    // std::cout<<"Matrix L final:";
    // print_matrix(L,nU,nF);
    // std::cout<<"Matrix R final:";
    // print_matrix(R,nF,nI);
    // std::cout<<"Matrix B final:";
    // print_matrix(B,nU,nI);

    double max = 0;
    double foo = 0;
    int col = 0;

    // Initilization file
    name_file = name_file+".out"; // add extension .out to final file
    myfile.open (name_file);
    myfile << "";
    myfile.close();

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

        // Write into file matFact inst0.in
        myfile.open (name_file, std::ios_base::app);
        myfile << col << "\n";
        myfile.close();
    }

  return 0;
}
