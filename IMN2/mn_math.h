#ifndef MNMATH_H
#define	MNMATH_H

#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>

using namespace std;


template <class T = double> 
class mnmath
{
public:



/* Alokuje tablice dwuwymiarowa o wymiarze [nrowa][ncols]
 * Indeksy tablicy zaczynaja sie od 0 */
static T** matrix(long nrh,  long nch)
{
    T** m;
    m = new T*[nrh];
    for(int i = 0 ; i < nrh ; i++) m[i] = new T[nch];
    for(int i = 0 ; i < nrh ; i++)
    for(int j = 0 ; j < nch ; j++) m[i][j] = 0;
	return m;
}
/* Ustawia wartosci macierzy na wartosc value*/
static void set_matrix(T **m,long nrh,  long nch,T value){
        for(int i = 0 ; i < nrh ; i++)
        for(int j = 0 ; j < nch ; j++) m[i][j] = value;
}

static void set_matrix_row(T **m,long row,long ncols,T value){
    for(int i = 0 ; i < ncols ; i++) {
        m[row][i] = value;
    }
}

/* Wypelnia wiersz row macierzy m wartoscia value */
static void set_matrix_row(T **m,long row,long ncols,T* values){
    for(int i = 0 ; i < ncols ; i++) {
        m[row][i] = values[i];
    }
}
/* Wypelnia kolumne col macierzy m wartoscia value */
static void set_matrix_col(T **m,long col,long nrows,T value){
    for(int i = 0 ; i < nrows ; i++) {
        m[i][col] = value;
    }
}

/* Wypelnia kolumne col macierzy m wartosciami values[i] */
static void set_matrix_col(T **m,long col,long nrows,T* values){
    for(int i = 0 ; i < nrows ; i++) {
        m[i][col] = values[i];
    }
}

/* Kopiuje wartosci macierzy srcM do dstM*/
static void copy_matrix(T **scrM, T **dstM,long nrh,long nch){
        for(int i = 0 ; i < nrh ; i++)
        for(int j = 0 ; j < nch ; j++) dstM[i][j] = scrM[i][j];
}

/* Wykonuje mnozenie macierzy razy wektor */
static void matrix_vector(T **mat, T *vec,T *vout,long nrh){
        for(int i = 0 ; i < nrh ; i++){
            vout[i] = 0;
            for(int j = 0 ; j < nrh ; j++)  vout[i] += mat[i][j] * vec[j];
        }
}

/* Wykonuje operacje mnozenia macierzy D =  A * B, przy czym wszystkie macierze sa kwadratowe */
static void matrix_matrix(T **matA, T **matB,T **matD,long nrh){
        for(int i = 0 ; i < nrh ; i++){
        for(int j = 0 ; j < nrh ; j++){
            matD[i][j] = 0;
            for(int k = 0 ; k < nrh ; k++){
                matD[i][j] += matA[i][k]*matB[k][j];
            }
        }}
}


/* Alokuje jednowymiarowy wektor o wymiarze [size] */
static T* vector(long nrh){
    T* v;
    v = new T[nrh];
    for(int i = 0 ; i < nrh ; i++) v[i] = 0;
    return v;
}

/* Ustawia wartosci wektora na value */
static void set_vector(T *v,long size,T value){
    for(int i = 0 ; i < size ; i++) v[i] = value;
}

/* Kopiuje wartosci wektora scrV do dstV */
static void copy_vector(T *scrV, T *dstV,long size){
    for(int i = 0 ; i < size ; i++) dstV[i] = scrV[i];
}

/* Zwalnia pamiec zajeta przez wektor*/
static void free_vector(T* v){
    delete [] v;
}
/* Zwalnia macierz utworzona za pomoca procedury matrix(..)*/
static void free_matrix(T **m,int nrows)
{
    for(int i = 0 ; i < nrows ; i++) delete [] m[i];
    delete [] m;
}


// -------------------------------------------------------------------------
//                      Funkcje wejscia / wyjscia
// -------------------------------------------------------------------------

/**
 * Zapisuje macierz mat do pliku o nazwie filename. Liczba kolumna i wierszy moze
 * byc mniejsza niz rzeczywista rozmiar macierzy.
 * @param forceEnter - dodaje pusta linie co forceEnter liczbe wierszy.
 */
static void write_matrix(const char* filename,T **mat,int nrows,int ncols,int forceEnter){
    FILE *f;
    int i,j;
    f = fopen(filename,"w");
    for(i = 0; i < nrows ; i++){
        for(j = 0; j < ncols ; j++){
			 
            if(sizeof(T)==sizeof(double) || sizeof(T)==sizeof(float)) fprintf(f,"%20.10e\t",double(mat[i][j]));
            //if(sizeof(T)==sizeof(int)) fprintf(f,"%16d\t",int(mat[i][j]));
        }
        fprintf(f,"\n");
        if(forceEnter!=0)if((i+1)%forceEnter == 0) fprintf(f,"\n");
    }
    fclose(f);
}


/**
 * Wypisuje macierz do standardowego wyjscia.
 * @param title - naglowek wypisywania
 * @param format - format zapisu w funkcji fprintf, np. dla double moze byc %16.6e.
 * @param mat   - wskaznik do macierzy
 * @param nrows - liczba wierszy
 * @param ncols - liczba kolumn
 */
static void print_matrix(const char* title,const char* format, T **mat,int nrows,int ncols){
    printf("%s\n",title);
    for(int i = 0; i < nrows ; i++){
        for(int j = 0; j < ncols ; j++){
            printf(format,double(mat[i][j]));
        }
        printf("\n");
    }
}


};


#endif	/* IMNMATH_H */

