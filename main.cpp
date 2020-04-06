/**********************************************************************************************************/
/** Programa que simula el modelo de Ising mediante el algoritmo de metropolis                           **/
/** El programa genera un archivo txt con la evoluc�n de la red de electrones a lo largo del tiempo para **/
/** una temperatura y un tama�o de red dados.                                                            **/
/**********************************************************************************************************/
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include "gsl_rng.h"

gsl_rng *tau;//Definici�n del puntero para generar n�meros aleatorios

using namespace std;

void inicializar (int N, double **M);
double exponencial (int m, int n, double T,int N,double **M);
double p (int a);

int main()
{
    int i,j,k;//Contadores
    int n,m;//posiciones aleatorias en la red
    double num;//Variable auxiliar para almacenar un numero
    double prob;//Variable para almacenar la probabilidad
    double T;//Variable temperatura
    int N;//Dimensiones de la red de electrones
    int semilla;//Semilla inicial para la generaci�n de n�meros aleatorios
    double **M;//Matriz en la que se almacena la red
    extern gsl_rng *tau;//Puntero de generaci�n de n�meros aleatorios
    ofstream fich;//Fichero para almacenar los datos



    semilla=1869732;//Semilla para la generaci�n de numeros aleatorios

    fich.open("Red semilla: "+to_string(semilla));

    //Inicializaci�n del generador de n�meros aleatorios
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    //Dimensiones de la red
    cout<<"Introduzca el numero de electrones en la red"<<endl;
    cin>>N;

    //Pido el vlor de la temperatura
    cout<<"Introduzca una temperatura"<<endl;
    cin>>T;

    //Definici�n de la memoria dinamica de la matriz
    M=new double*[N];
    for(i=0;i<N;i++)
        M[i]=new double[N];

    //Inicializaci�n de los valores para el spin de la red con condiciones de contorno periodicas
    inicializar(N,M);

    //Almaceno la primera configuraci�n en un fichero
    for(i=0;i<N;i++)
    {
        for(k=0;k<N;k++)
        {
            fich<<M[i][k]<<"\t";
        }
        fich<<endl;
    }
    fich<<endl<<endl;

    //Implementaci�n del algoritmo de montecarlo para simular el modelo de Ising
    //El algoritmo dar�  pasos montecarlo cambiando los valores de los espines
    for (j=0;j<500;j++)
    {

        for (i=0;i<(N*N);i++)
        {
            //Generaci�n de posiciones aleatorias (paso 1)
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N);

            //Evaluaci�n de la probabilidad (paso 2)
            num=exponencial(m,n,T,N,M);
            prob=p(num);

            //Generaci�n de un n�mero al azar y comprobaci�n de si queda dentro de la probabilidad (paso 3)
            num=gsl_rng_uniform(tau);
            if (num<prob)
                M[m][n]=-1*M[m][n];
        }

        //Almaceno los resultados en un fichero tras cada paso montecarlo
        for(i=0;i<N;i++)
        {
            for(k=0;k<N;k++)
            {
                fich<<M[i][k]<<"\t";
            }
            fich<<endl;
        }

        fich<<endl<<endl;

    }


    //Borrado de la memoria dinamica
    for (i=0;i<N;i++)
        delete[] M[i];
    delete[] M;

    //Cierre de ficheros
    fich.close();

    return 0;
}


//Funci�n que inicializa la red con condiciones de contorno periodicas
void inicializar (int N, double **M)
{
    int i,j;
    double num; //Variable auxiliar para almacenar n�meros
    extern gsl_rng *tau;//Puntero de generaci�n de n�meros aleatorios

    //Inicializaci�n de las condiciones de contorno
    for(j=0;j<N;j++)
    {
        num=gsl_rng_uniform(tau);
        if (num>0.5)
        {
            M[0][j]=1;
            M[N-1][j]=1;
        }

        else
        {
            M[0][j]=-1;
            M[N-1][j]=-1;
        }
    }

    for(i=1;i<N;i++)
    {
        num=gsl_rng_uniform(tau);
        if (num>0.5)
        {
            M[i][0]=1;
            M[i][N-1]=1;
        }

        else
        {
            M[i][0]=-1;
            M[i][N-1]=-1;
        }
    }

    //inicializaci�n del resto de la red
    for (i=1;i<(N-1);i++)
        for (j=1;j<(N-1);j++)
        {
            num=gsl_rng_uniform(tau);
            if (num>0.5)
                M[i][j]=1;
            else
                M[i][j]=-1;
        }

    return;
}

//Funci�n que calcula la probabilidad de encontrar una configuraci�n espines en la red
//La funcion tiene en cuenta las condiciones de contorno periodicas para calcular la probabilidad
double exponencial (int m, int n, double T,int N,double **M)
{
    if ((n-1)==-1)
    {
        if ((m-1)==-1)
            return exp(-(2*M[m][n]*(M[m+1][n]+M[0][n]+M[m][n+1]+M[m][0]))/T);

        else
        {
            if ((m+1)==N)
                return exp(-(2*M[m][n]*(M[0][n]+M[m-1][n]+M[m][n+1]+M[m][0]))/T);

            else
                return exp(-(2*M[m][n]*(M[m+1][n]+M[m-1][n]+M[m][n+1]+M[m][0]))/T);
        }
    }




    else
    {
        if ((m-1)==-1)
        {
            if ((n+1)==N)
                return exp(-(2*M[m][n]*(M[m+1][n]+M[0][n]+M[m][0]+M[m][n-1]))/T);

            else
                return exp(-(2*M[m][n]*(M[m+1][n]+M[0][n]+M[m][n+1]+M[m][n-1]))/T);
        }
        else
        {
            if ((n+1)==N)
            {
                if ((m+1)==N)
                    return exp(-(2*M[m][n]*(M[0][n]+M[m-1][n]+M[m][0]+M[m][n-1]))/T);

                else
                    return exp(-(2*M[m][n]*(M[m+1][n]+M[m-1][n]+M[m][0]+M[m][n-1]))/T);
            }

            else
            {
                if ((m+1)==N)
                    return exp(-(2*M[m][n]*(M[0][n]+M[m-1][n]+M[m][n+1]+M[m][n-1]))/T);

                else
                    return exp(-(2*M[m][n]*(M[m+1][n]+M[m-1][n]+M[m][n+1]+M[m][n-1]))/T);
            }
        }
    }
}

//Funcion que evalua si la probabilidad calculada por la funci�n anterior es suficientemente posible mediante el uso de otro numero generado aleatoriamente
double p (int a)
{
    if (a<1)
        return a;
    else
        return 1;
}
