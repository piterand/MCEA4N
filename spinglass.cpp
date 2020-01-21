#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <mpi.h>


// NO PPF!!!!!!!!!!
////////////////////
#define n 100		 	///Ввод размерности решетки
#define prohod_MC2 n*n		//количество проходов внутри
float Polevar = 0;           // Поле
unsigned long long  Prohod_MC_equilibration = 1e5;		/// Ввод числа шагов разогрева
unsigned long long  Prohod_MC_sampling = 1e6;		/// Ввод числа шагов Монте-Карло
unsigned long long between_exchange_equilibration=1e3; //между репличными обменами
unsigned long long between_exchange=1e2; //между репличными обменами
bool replica_on=1;
int current_value = 0;              // выводить текущие значения или нет
unsigned long long CountOut = 1e5;                             // количество вывода

///////////////////////

using namespace std;

///////////////////////////////////		Ввод переменных
int choice_linear_dimension;			// выбор линейного размера системы
double choice_step_temperature;			// выбор шага температуры
int num_podhod=1000;			    		// количество попыток найти систему связей с мин. кол-м ошибок.
double mintemp=3.1;						// выбор min температуры
double maxtemp=4.0;						// выбор max температуры



////////////////////////////////////
void Randsvaz(short *svaz);
void Vivoderrors(int num, short *svaz);
void Slspinn (short *spin);
void SpinFerr(short *spin);
void Vivodspin(int spin);
int  Slsvaz( short *svaz, int num_podhod);
void VivodSSS( short *spin,short *svaz);
void Energy(short *spin,short *svaz,float *E);
float EnergyIJ(int i,int j,short *spin,short *svaz);
void VivodEnergy( float *E);
void CopyE(float *E,float *E1);
int Claster( float *E1);
int ClasterFerr( float *E1);
void ColorOut2( float *E, short *spin);
void ColorOut( float *E, short *spin);
float SumEnergy(float *E);
void ColorOutIJ( float *E, short *spin, int i,int j);
void ColorOutSpin( float *E, short *spin);
int MaxClass(int per, float *E1, int tmp, int *Ochered);
int MaxClassFerr(int per, float *E1, int tmp, int *Ochered);
int sosed_sv(int top);
int sosed_sn(int top);
int sosed_sl(int top);
int sosed_sp(int top);
double SROTKL(double OutputMass, double SRKVOTKL);
double MagnMet( short *spin);
void PhazDiagramm(int var04, int var02, int var00, int var20, int var40, float *E);
void MCL(double mintemp, double maxtemp, short *spin, short *svaz, float *E, float *E1, int NUM_step_mcl);
void MCProhod(double temp,short *spin,short *svaz,float *E);

////////////////////////////////////////////////////////вывод координат ошибок
void Vivoderrors(int num, short *svaz)
        {
            int x = num;
            cout<<"Координаты ошибок в решетке:"<<endl;
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < x; j++)
                {
                    if ((svaz[i, j] + svaz[i, j + 1] + svaz[i + num, j] + svaz[i + num + 1, j]) != 0)//////////// вывод ошибок
                    {
                        cout<<svaz[i*(num+1) + j]<<",";
                        cout<<svaz[(i+num)*(num+1) + j]<<endl;
                        cout<<svaz[i*(num+1) + (j+1)]<<",";
                        cout<<svaz[(i+num+1)*(num+1) + j]<<endl;                        
                    }

                }
            }
            cout<<endl;
        }

///////////////////////////////////функция визуализации спинов
// static void VivodSSS(int num, short *spin, short *svaz)
//         {
//             cout<<"Вывод спинов со связями"<<end;
//             int x = num;
//             for (int i = 0; i < x; i++)
//             {
//                 for (int j = 0; j < (2 * x); j++)
//                 {
//                     if (j < x)
//                     {
//                         if (svaz[(i)*(num+1)] == 1) //if (svaz[x + i, j] == 1)
//                             cout<<" +";
//                         else
//                             cout<<" -";
//                     }
//                     if (j == x)
//                         cout<endl;
//                     if (j >= x)
//                     {
//                         if (svaz[i, j - x] == 1) //svaz[i, j - x] == 1
//                             cout<<" +";
//                         else
//                             cout<<" -";

//                         if (spin[i, j - x] == -1) //spin[i, j - x] == -1
//                         {
//                             cout<<"↓";
//                         }
//                         else
//                         {
//                             cout<<"↑";
//                         }
//                         if (j == (2 * x - 1)) //j == (2 * x - 1)
//                         {
//                             if (svaz[i, x] == 1)
//                                 cout<<" +";
//                             else
//                                 cout<<" -";
//                         }
//                     }
//                 }
//                 count<<endl;
//             }
//             for (int j = 0; j < x; j++)
//             {
//                 if (svaz[x, j] == 1)
//                     cout<<" +";
//                 else
//                     cout<<" -";
//             }
//         }
void VivodSSS(int num, short *spin, short *svaz){
    cout<<"Вывод спинов со связями"<<endl;
        int x = num;
        int i,j;
        int koeff=0;

        for (i = 0; i < 2*x+1; i++)
        {
            if(!(i%2))
            {
                for(j=0; j<x; j++)
                {
                    if(!((x*i/2+j+1+koeff) % (x+1)))
                    	koeff++;

                    cout<<" ";
                    if(svaz[x*(x+(i+2)/2)+j+koeff]==1)
                        cout<<"+";
                    else
                        cout<<"-";
                }
                cout<<endl;
            }
            else
            {
                for(j=0; j<x; j++)
                {
                    if(svaz[(x+1)*(i-1)/2+j]==1)
                        cout<<"+";
                    else
                        cout<<"-";

                    if(spin[x*(i-1)/2+j]==1)
                        cout<<"↑";
                    else
                        cout<<"↓";
                }

                if(svaz[(x+1)*(i-1)/2+j]==1)
                    cout<<"+";
                else
                    cout<<"-";

                cout<<endl;
            }
        }
        // for(i=0; i<(x*2+1)*(x+1); i++)
        //     cout<<svaz[i];
}
//////////////////////////////////////////////////// Забивает матрицу спинов 1

void SpinFerr(short *spin)	
{
    for(int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            spin[i*n+j]=1;
        }
    }
}

void Slspinn(short *spin)
{
    for(int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            spin[i*n+j]=(rand()%2)*2-1;

        }
    }
}

void Randsvaz(short *svaz)
{
    for(int i=0;i<(n*2+1)*(n+1);i++)
    {
        svaz[i]=(rand()%2)*2-1;
    }
}

int Slsvaz(short *svaz,int num_podhod)
{
    int x = n;
    int errorN=1;
    int min_error =1000;
    int l=0;
    short *svaz_min_error = new short[(n*2+1)*(n+1)];
    ////////////////////////////////////////////////
    while(errorN != 0 && l < num_podhod)
    {
        if (l!=0)
        {
            if (errorN < min_error)
            {
                min_error = errorN;
                for (int i = 0; i < (n*2+1)*(n+1); i++)
                    svaz_min_error[i] = svaz[i];
                                    
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    svaz[i*(n+1) + j] = 0;
                    svaz[(i+n)*(n+1) + j] = 0;
                    svaz[i*(n+1) + (j+1)] = 0;
                    svaz[(i+n+1)*(n+1) + j] = 0;
                }
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    svaz[i*(n+1) + j] = 0;
                    svaz[(i+n)*(n+1) + j] = 0;
                    svaz[i*(n+1) + (j+1)] = 0;
                    svaz[(i+n+1)*(n+1) + j] = 0;
                }
            }
        }
        errorN = 0;

        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < x; j++)   ////////// Забиваем массив со связями сл.величинами так, что бы в сумме около одного спина был 0
            {
                int Sl1;
                int per = 0; ;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 0;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 1;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 2;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 3;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 4;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 5;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 6;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 7;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 8;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 9;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 10;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0)
                    per = 11;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 12;
                if (svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 13;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 14;
                if (svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0)
                    per = 15;
                
                    switch (per)
                {
                    case 0:                 //////////  svaz[i][j] == 0 && svaz[i][j+1] == 0 && svaz[i+n][j] == 0 && svaz[i+n+1][j] == 0

                        if (j != x - 2)
                        {
                            svaz[(i)*(n+1) + (j)] = rand()%(2);
                            if (svaz[(i)*(n+1) + (j)] == 0)
                                svaz[(i)*(n+1) + (j)] = -1;
                            Sl1 = rand()%(3);

                            switch (Sl1)
                            {
                            case 0:
                                svaz[(i)*(n+1) + (j+1)] = svaz[(i)*(n+1) + (j)];
                                if (svaz[(i)*(n+1) + (j+1)] == 1)
                                {
                                    svaz[(i+x)*(n+1) + (j)] = -1;
                                    svaz[(i+x+1)*(n+1) + (j)] = -1;
                                }
                                else
                                {
                                    svaz[(i+x)*(n+1) + (j)] = 1;
                                    svaz[(i+x+1)*(n+1) + (j)] = 1;
                                }
                                break;

                            case 1:
                                svaz[(i+x)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                if (svaz[(i+x)*(n+1) + (j)] == 1)
                                {
                                    svaz[(i)*(n+1) + (j+1)] = -1;
                                    svaz[(i+x+1)*(n+1) + (j)] = -1;
                                }
                                else
                                {
                                    svaz[(i)*(n+1) + (j+1)] = 1;
                                    svaz[(i+x+1)*(n+1) + (j)] = 1;
                                }
                                break;

                            case 2:
                                svaz[(i+x+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                if (svaz[(i+x+1)*(n+1) + (j)] == 1)
                                {
                                    svaz[(i)*(n+1) + (j+1)] = -1;
                                    svaz[(i+x)*(n+1) + (j)] = -1;
                                }
                                else
                                {
                                    svaz[(i)*(n+1) + (j+1)] = 1;
                                    svaz[(i+x)*(n+1) + (j)] = 1;
                                }

                                break;

                            }
                        }

                        if (j == x - 2)
                        {

                        }

                        

                        break;





                            case 1:               ////////////// svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0
                                Sl1 = rand()%(3);

                                switch (Sl1)
                                {
                                    case 0:
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i)*(n+1) + (j)];
                                        if (svaz[(i)*(n+1) + (j+1)] == 1)
                                        {
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                            svaz[(i+x+1)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                            svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 1:
                                        svaz[(i+x)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                        if (svaz[(i+x)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                            svaz[(i+x+1)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                            svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 2:
                                        svaz[(i+x+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                        if (svaz[(i+x+1)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                        }
                                        break;
                                }

                                break;



                            case 2:             ////////////// svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0
                                Sl1 = rand()%(3);

                                switch (Sl1)
                                {
                                    case 0:
                                        svaz[(i)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        if (svaz[(i)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                            svaz[(i+x+1)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                            svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 1:
                                        svaz[(i+x)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        if (svaz[(i+x)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j)] = -1;
                                            svaz[(i+x+1)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j)] = 1;
                                            svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 2:
                                        svaz[(i+x+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        if (svaz[(i+x+1)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j)] = -1;
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j)] = 1;
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                        }
                                        break;
                                }

                                break;

                            case 3:             //////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0

                                if (i == x - 2)
                                {
                                    svaz[(i)*(n+1) + (0)] = -svaz[(i + x)*(n+1) + (x - 1)];
                                    svaz[(i + 1)*(n+1) + (0)] = -svaz[(i + x + 2)*(n+1) + (x - 1)];
                                    if (svaz[(i)*(n+1) + ( 0)] + svaz[(i + x)*(n+1) + ( 0)] == 2 && svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] == -2 || svaz[(i)*(n+1) + ( 0)] + svaz[(i + x)*(n+1) + ( 0)] == -2 && svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] == 2)
                                    {
                                        svaz[(i + 1)*(n+1) + ( 0)] = -svaz[(i + x + 2)*(n+1) + ( 0)];
                                    }

                                    if (svaz[(i + n)*(n+1) + ( 0)] == svaz[(i + n + 2)*(n+1) + ( 0)])
                                        svaz[(i + x + 1)*(n+1) + ( 0)] = -svaz[(i + x)*(n+1) + ( 0)];
                                    else
                                    {
                                        if (svaz[(i)*(n+1) + ( 0)] + svaz[(i + n)*(n+1) + ( 0)] == 2 || svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] == 2)
                                            svaz[(i + x + 1)*(n+1) + ( 0)] = -1;
                                        if (svaz[(i)*(n+1) + ( 0)] + svaz[(i + n)*(n+1) + ( 0)] == -2 || svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] == -2)
                                            svaz[(i + x + 1)*(n+1) + ( 0)] = 1;
                                        if (svaz[(i)*(n+1) + ( 0)] + svaz[(i + n)*(n+1) + ( 0)] == 0 && svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] == 0)
                                            svaz[(i + x + 1)*(n+1) + ( 0)] = 1;
                                    }
                                    if (svaz[(i)*(n+1) + ( 0)] + svaz[(i + x)*(n+1) + ( 0)] + svaz[(i + x + 1)*(n+1) + ( 0)] > 0)
                                        svaz[(i)*(n+1) + ( 1)] = -1;
                                    else
                                        svaz[(i)*(n+1) + ( 1)] = 1;

                                    if (svaz[(i + 1)*(n+1) + ( 0)] + svaz[(i + x + 1)*(n+1) + ( 0)] + svaz[(i + x + 2)*(n+1) + ( 0)] > 0)
                                        svaz[(i + 1)*(n+1) + ( 1)] = -1;
                                    else
                                        svaz[(i + 1)*(n+1) + ( 1)] = 1;



                                    for (int r = 1; r < x - 1; r++)
                                    {
                                        if (svaz[(i)*(n+1) + ( r)] + svaz[(i + x)*(n+1) + ( r)] == 2 || svaz[(i + 1)*(n+1) + ( r)] + svaz[(i + x + 2)*(n+1) + ( r)] == 2)
                                            svaz[(i + x + 1)*(n+1) + ( r)] = -1;
                                        else
                                            svaz[(i + x + 1)*(n+1) + ( r)] = 1;

                                        if (svaz[(i)*(n+1) + ( r)] + svaz[(i + x)*(n+1) + ( r)] + svaz[(i + x + 1)*(n+1) + ( r)] > 0)
                                            svaz[(i)*(n+1) + ( r + 1)] = -1;
                                        else
                                            svaz[(i)*(n+1) + ( r + 1)] = 1;

                                        if (svaz[(i + 1)*(n+1) + ( r)] + svaz[(i + x + 1)*(n+1) + ( r)] + svaz[(i + x + 2)*(n+1) + ( r)] > 0)
                                            svaz[(i + 1)*(n+1) + ( r + 1)] = -1;
                                        else
                                            svaz[(i + 1)*(n+1) + ( r + 1)] = 1;
                                        int w = r;
                                        /*
                                        if (svaz[(i)*(n+1) + ( r)] + svaz[(i + x)*(n+1) + ( r)] == 2 && svaz[(i + 1)*(n+1) + ( r)] + svaz[(i + x + 2)*(n+1) + ( r)] == -2 || svaz[(i)*(n+1) + ( r)] + svaz[(i + x)*(n+1) + ( r)] == -2 && svaz[(i + 1)*(n+1) + ( r)] + svaz[(i + x + 2)*(n+1) + ( r)] == 2)

                                        {
                                            while (svaz[(i)*(n+1) + ( w)] + svaz[(i + x)*(n+1) + ( w)] == 0 && svaz[(i + 1)*(n+1) + ( w)] + svaz[(i + x + 2)*(n+1) + ( w)] == 0)
                                            {
                                                w--;
                                            }
                                            svaz[(i + x + 1)*(n+1) + ( w)] = -svaz[(i + x + 1)*(n+1) + ( w)];
                                            svaz[(i)*(n+1) + ( w + 1)] = -svaz[(i)*(n+1) + ( w + 1)];
                                            svaz[(i + 1)*(n+1) + ( w + 1)] = svaz[(i + 1)*(n+1) + ( w + 1)];
                                            w++;
                                            for (; w <= r; w++)
                                            {
                                                
                                                
                                                if (svaz[(i)*(n+1) + ( w)] + svaz[(i + x)*(n+1) + ( w)] == 2 || svaz[(i + 1)*(n+1) + ( w)] + svaz[(i + x + 2)*(n+1) + ( w)] == 2)
                                                    svaz[(i + x + 1)*(n+1) + ( w)] = -1;
                                                else
                                                    svaz[(i + x + 1)*(n+1) + ( w)] = 1;

                                                if (svaz[(i)*(n+1) + ( w)] + svaz[(i + x)*(n+1) + ( w)] + svaz[(i + x + 1)*(n+1) + ( w)] > 0)
                                                    svaz[(i)*(n+1) + ( w + 1)] = -1;
                                                else
                                                    svaz[(i)*(n+1) + ( w + 1)] = 1;

                                                if (svaz[(i + 1)*(n+1) + ( w)] + svaz[(i + x + 1)*(n+1) + ( w)] + svaz[(i + x + 2)*(n+1) + ( w)] > 0)
                                                    svaz[(i + 1)*(n+1) + ( w + 1)] = -1;
                                                else
                                                    svaz[(i + 1)*(n+1) + ( w + 1)] = 1;
                                            }
                                        }
                                        */





                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(3);

                                    switch (Sl1)
                                    {
                                        case 0:
                                            svaz[(i)*(n+1) + (j)] = svaz[(i+x)*(n+1) + (j)];
                                            if (svaz[(i)*(n+1) + (j)] == 1)
                                            {
                                                svaz[(i)*(n+1) + (j+1)] = -1;
                                                svaz[(i+x+1)*(n+1) + (j)] = -1;
                                            }
                                            else
                                            {
                                                svaz[(i)*(n+1) + (j+1)] = 1;
                                                svaz[(i+x+1)*(n+1) + (j)] = 1;
                                            }
                                            break;

                                        case 1:
                                            svaz[(i)*(n+1) + (j+1)] = svaz[(i+x)*(n+1) + (j)];
                                            if (svaz[(i)*(n+1) + (j+1)] == 1)
                                            {
                                                svaz[(i)*(n+1) + (j)] = -1;
                                                svaz[(i+x+1)*(n+1) + (j)] = -1;
                                            }
                                            else
                                            {
                                                svaz[(i)*(n+1) + (j)] = 1;
                                                svaz[(i+x+1)*(n+1) + (j)] = 1;
                                            }
                                            break;

                                        case 2:
                                            svaz[(i+x+1)*(n+1) + (j)] = svaz[(i+x)*(n+1) + (j)];
                                            if (svaz[(i+x+1)*(n+1) + (j)] == 1)
                                            {
                                                svaz[(i)*(n+1) + (j)] = -1;
                                                svaz[(i)*(n+1) + (j+1)] = -1;
                                            }
                                            else
                                            {
                                                svaz[(i)*(n+1) + (j)] = 1;
                                                svaz[(i)*(n+1) + (j+1)] = 1;
                                            }
                                            break;
                                    }
                                }

                                break;

                            case 4:    ////////////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0


                                Sl1 = rand()%(3);

                                switch (Sl1)
                                {
                                    case 0:
                                        svaz[(i)*(n+1) + (j)] = svaz[(i+x+1)*(n+1) + (j)];
                                        if (svaz[(i)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 1:
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i+x+1)*(n+1) + (j)];
                                        if (svaz[(i)*(n+1) + (j+1)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j)] = -1;
                                            svaz[(i+x)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j)] = 1;
                                            svaz[(i+x)*(n+1) + (j)] = 1;
                                        }
                                        break;

                                    case 2:
                                        svaz[(i+x)*(n+1) + (j)] = svaz[(i+x+1)*(n+1) + (j)];
                                        if (svaz[(i+x)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j)] = -1;
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j)] = 1;
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                        }
                                        break;
                                }

                                break;



                            case 5:         ////////////////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] == 0


                                if (svaz[(i)*(n+1) + (j)] == svaz[(i)*(n+1) + (j+1)])
                                {
                                    if (svaz[(i)*(n+1) + (j)] == 1)
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = -1;
                                        svaz[(i+n+1)*(n+1) + (j)] = -1;
                                    }
                                    else
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = 1;
                                        svaz[(i+n+1)*(n+1) + (j)] = 1;
                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(2);
                                    if (Sl1 == 0)
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                        svaz[(i+n+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                    }
                                    else
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        svaz[(i+n+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                    }

                                }

                                break;

                            case 6:         ////////////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0

                                if (j == x - 2 && (svaz[(i)*(n+1) + ( j + 2)] + svaz[(i + x)*(n+1) + ( j + 1)] == 2 || svaz[(i)*(n+1) + ( j + 2)] + svaz[(i + x)*(n+1) + ( j + 1)] == -2))
                                {
                                    if (svaz[(i)*(n+1) + (j)] + svaz[(i+n)*(n+1) + (j)] != -(svaz[(i)*(n+1) + ( j + 2)] + svaz[(i + x)*(n+1) + ( j + 1)]))
                                    {
                                        if (svaz[(i)*(n+1) + ( j + 2)] + svaz[(i + x)*(n+1) + ( j + 1)] == 2)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                            if (svaz[(i)*(n+1) + (j)] + svaz[(i+x)*(n+1) + (j)] + svaz[(i)*(n+1) + (j+1)] > 0)
                                                svaz[(i+x+1)*(n+1) + (j)] = -1;
                                            else
                                                svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                        if (svaz[(i)*(n+1) + ( j + 2)] + svaz[(i + x)*(n+1) + ( j + 1)] == -2)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                            if (svaz[(i)*(n+1) + (j)] + svaz[(i+x)*(n+1) + (j)] + svaz[(i)*(n+1) + (j+1)] > 0)
                                                svaz[(i+x+1)*(n+1) + (j)] = -1;
                                            else
                                                svaz[(i+x+1)*(n+1) + (j)] = 1;
                                        }
                                    }
                                    else
                                    {
                                        int p = j;

                                        while (p >= 1 && svaz[(i)*(n+1) + ( p)] + svaz[(i + x)*(n+1) + ( p)] != 0)
                                            p--;

                                        p++;
                                        for (; p < x; p++)
                                        {
                                            svaz[(i)*(n+1) + ( p)] = -svaz[(i + x)*(n+1) + ( p)];
                                            svaz[(i + x + 1)*(n+1) + ( p - 1)] = -svaz[(i)*(n+1) + ( p)];
                                            if (p == x - 1)
                                                svaz[(i + x + 1)*(n+1) + ( p)] = -svaz[(i)*(n+1) + ( p + 1)];
                                        }
                                    }

                                }
                                else
                                {

                                    if (svaz[(i)*(n+1) + (j)] == svaz[(i+n)*(n+1) + (j)])
                                    {
                                        if (svaz[(i)*(n+1) + (j)] == 1)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = -1;
                                            svaz[(i+n+1)*(n+1) + (j)] = -1;
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = 1;
                                            svaz[(i+n+1)*(n+1) + (j)] = 1;
                                        }
                                    }
                                    else
                                    {
                                        Sl1 = rand()%(2);
                                        if (Sl1 == 0)
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = svaz[(i)*(n+1) + (j)];
                                            svaz[(i+n+1)*(n+1) + (j)] = svaz[(i+n)*(n+1) + (j)];
                                        }
                                        else
                                        {
                                            svaz[(i)*(n+1) + (j+1)] = svaz[(i+n)*(n+1) + (j)];
                                            svaz[(i+n+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                        }

                                    }

                                }


                                break;

                            case 7:         ////////////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0

                                if (svaz[(i)*(n+1) + (j)] == svaz[(i+n+1)*(n+1) + (j)])
                                {
                                    if (svaz[(i)*(n+1) + (j)] == 1)
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = -1;
                                        svaz[(i)*(n+1) + (j+1)] = -1;
                                    }
                                    else
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = 1;
                                        svaz[(i)*(n+1) + (j+1)] = 1;
                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(2);
                                    if (Sl1 == 0)
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i)*(n+1) + (j)];
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i+n+1)*(n+1) + (j)];
                                    }
                                    else
                                    {
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i+n+1)*(n+1) + (j)];
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i)*(n+1) + (j)];
                                    }

                                }

                                break;


                            case 8:             /////////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0

                                if (svaz[(i)*(n+1) + (j+1)] == svaz[(i+n)*(n+1) + (j)])
                                {
                                    if (svaz[(i)*(n+1) + (j+1)] == 1)
                                    {
                                        svaz[(i)*(n+1) + (j)] = -1;
                                        svaz[(i+n+1)*(n+1) + (j)] = -1;
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = 1;
                                        svaz[(i+n+1)*(n+1) + (j)] = 1;
                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(2);
                                    if (Sl1 == 0)
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i+n)*(n+1) + (j)];
                                        svaz[(i+n+1)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        svaz[(i+n+1)*(n+1) + (j)] = svaz[(i+n)*(n+1) + (j)];
                                    }

                                }

                                break;

                            case 9:             /////////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0

                                if (svaz[(i)*(n+1) + (j+1)] == svaz[(i+n+1)*(n+1) + (j)])
                                {
                                    if (svaz[(i)*(n+1) + (j+1)] == 1)
                                    {
                                        svaz[(i)*(n+1) + (j)] = -1;
                                        svaz[(i+n)*(n+1) + (j)] = -1;
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = 1;
                                        svaz[(i+n)*(n+1) + (j)] = 1;
                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(2);
                                    if (Sl1 == 0)
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i+n+1)*(n+1) + (j)];
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i)*(n+1) + (j+1)];
                                        svaz[(i+n)*(n+1) + (j)] = svaz[(i+n+1)*(n+1) + (j)];
                                    }

                                }

                                break;

                            case 10:             /////////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0

                                if (svaz[(i+x)*(n+1) + (j)] == svaz[(i+n+1)*(n+1) + (j)])
                                {
                                    if (svaz[(i+x)*(n+1) + (j)] == 1)
                                    {
                                        svaz[(i)*(n+1) + (j)] = -1;
                                        svaz[(i)*(n+1) + (j+1)] = -1;
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = 1;
                                        svaz[(i)*(n+1) + (j+1)] = 1;
                                    }
                                }
                                else
                                {
                                    Sl1 = rand()%(2);
                                    if (Sl1 == 0)
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i+n+1)*(n+1) + (j)];
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i+x)*(n+1) + (j)];
                                    }
                                    else
                                    {
                                        svaz[(i)*(n+1) + (j)] = svaz[(i+x)*(n+1) + (j)];
                                        svaz[(i)*(n+1) + (j+1)] = svaz[(i+n+1)*(n+1) + (j)];
                                    }
                                }

                                break;


                            case 11:  ////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] == 0
                                if (svaz[(i)*(n+1) + (j)] + svaz[(i)*(n+1) + (j+1)] + svaz[(i+n)*(n+1) + (j)] > 0)
                                    svaz[(i+n+1)*(n+1) + (j)] = -1;
                                else
                                    svaz[(i+n+1)*(n+1) + (j)] = 1;


                                break;

                            case 12:     ////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] == 0 && svaz[(i+n+1)*(n+1) + (j)] != 0
                                if (svaz[(i)*(n+1) + (j)] + svaz[(i)*(n+1) + (j+1)] + svaz[(i+n+1)*(n+1) + (j)] > 0)
                                    svaz[(i+n)*(n+1) + (j)] = -1;
                                else
                                    svaz[(i+n)*(n+1) + (j)] = 1;

                                break;


                            case 13:     ////////svaz[(i)*(n+1) + (j)] == 0 && svaz[(i)*(n+1) + (j+1)] != 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0
                                if (svaz[(i)*(n+1) + (j+1)] + svaz[(i+n)*(n+1) + (j)] + svaz[(i+n+1)*(n+1) + (j)] > 0)
                                    svaz[(i)*(n+1) + (j)] = -1;
                                else
                                    svaz[(i)*(n+1) + (j)] = 1;

                                break;

                            case 14:     ////////svaz[(i)*(n+1) + (j)] != 0 && svaz[(i)*(n+1) + (j+1)] == 0 && svaz[(i+n)*(n+1) + (j)] != 0 && svaz[(i+n+1)*(n+1) + (j)] != 0
                                if (svaz[(i)*(n+1) + (j)] + svaz[(i+n)*(n+1) + (j)] + svaz[(i+n+1)*(n+1) + (j)] > 0)
                                    svaz[(i)*(n+1) + (j+1)] = -1;
                                else
                                    svaz[(i)*(n+1) + (j+1)] = 1;

                                break;

                            case 15:
                                break;
                        }

                                if (j == 0)  ////////////////// если спин крайний, то выполнение граничных условий:
                                {
                                    svaz[(i)*(n+1) + (x)] = svaz[(i)*(n+1) + (0)];

                                }
                                if (i == 0)
                                {
                                    svaz[(x * 2)*(n+1) + (j)] = svaz[(x)*(n+1) + (j)];
                                }


                        

                    }
                }


                for (int i = 0; i < x; i++)
                {
                    for (int j = 0; j < x; j++)
                    {
                        if (svaz[(i)*(n+1) + (j)] + svaz[(i)*(n+1) + (j + 1)] + svaz[(i + n)*(n+1) + (j)] + svaz[(i + n + 1)*(n+1) + (j)] != 0)//////////// подсчет  ошибок
                        {
                            errorN++;

                        }
                        
                    }
                }

                l++;
            }
            if (errorN != 0)
            {
                for (int i = 0; i < 2 * x + 1; i++)
                {
                    for (int j = 0; j < x + 1; j++)
                    {
                        svaz[(i)*(n+1) + (j)] = svaz_min_error[(i)*(n+1) + (j)];
                    }
                }
            }
            else min_error = 0;
            delete[] svaz_min_error;
            return min_error;
}
///////////////////////////////////////////////////////Подсчитывает энегрию системы
void Energy(short *spin,short *svaz,float *E)
{
    int x = n;
    E[0*n+0] = -svaz[(0)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(0)*n + ( x - 1)] - svaz[(0)*(n+1) + ( 1)] * spin[(0)*n + ( 0)] * spin[(0)*n + ( 1)] - svaz[(x)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(x - 1)*n + ( 0)] - svaz[(x + 1)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(1)*n + ( 0)] - Polevar;
    E[0*n+(x - 1)] = -svaz[(0)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(0)*n + ( x - 2)] - svaz[(0)*(n+1) + ( x)] * spin[(0)*n + ( x - 1)] * spin[(0)*n + ( 0)] - svaz[(x)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] - svaz[(x + 1)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(1)*n + ( x - 1)] - Polevar;
    E[(x - 1)*n + 0] = -svaz[(x - 1)*(n+1) + ( 0)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( 0)] - svaz[(x - 1)*(n+1) + ( 1)] * spin[(x - 1)*n + ( 0)] * spin[(x - 1)*n + ( 1)] - svaz[(2 * x - 1)*(n+1) + ( 0)] * spin[(x - 2)*n + ( 0)] * spin[(x - 1)*n + ( 0)] - svaz[(2 * x)*(n+1) + ( 0)] * spin[(x - 1)*n + ( 0)] * spin[(0)*n + ( 0)] - Polevar;
    E[(x - 1)*n + (x - 1)] = -svaz[(x - 1)*(n+1) + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 2)] - svaz[(x - 1)*(n+1) + ( x)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( 0)] - svaz[(2 * x - 1)*(n+1) + ( x - 1)] * spin[(x - 2)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] - svaz[(2 * x)*(n+1) + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] * spin[(0)*n + ( x - 1)] - Polevar;
    for (int i = 1; i < x - 1; i++)
    {
        for (int j = 1; j < x - 1; j++)
            E[i*n + j] = -svaz[(i)*(n+1) + ( j)] * spin[(i)*n + ( j - 1)] * spin[(i)*n + ( j)] - svaz[(i)*(n+1) + ( j + 1)] * spin[(i)*n + ( j)] * spin[(i)*n + ( j + 1)] - svaz[(i + x)*(n+1) + ( j)] * spin[(i - 1)*n + ( j)] * spin[(i)*n + ( j)] - svaz[(i + x + 1)*(n+1) + ( j)] * spin[(i)*n + ( j)] * spin[(i + 1)*n + ( j)] - Polevar;
    }
    for (int j = 1; j < x - 1; j++)
        E[0*n + j] = -svaz[(0)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(0)*n + ( j - 1)] - svaz[(0)*(n+1) + ( j + 1)] * spin[(0)*n + ( j)] * spin[(0)*n + ( j + 1)] - svaz[(x)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(x - 1)*n + ( j)] - svaz[(x + 1)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(1)*n + ( j)] - Polevar;
    for (int i = 1; i < x - 1; i++)
        E[i*n + 0] = -svaz[(i)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i)*n + ( x - 1)] - svaz[(i)*(n+1) + ( 1)] * spin[(i)*n + ( 0)] * spin[(i)*n + ( 1)] - svaz[(i + x)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i - 1)*n + ( 0)] - svaz[(i + x + 1)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i + 1)*n + ( 0)] - Polevar;
    for (int j = 1; j < x - 1; j++)
        E[(x - 1)*n + j] = -svaz[(x - 1)*(n+1) + ( j)] * spin[(x - 1)*n + ( j - 1)] * spin[(x - 1)*n + ( j)] - svaz[(x - 1)*(n+1) + ( j + 1)] * spin[(x - 1)*n + ( j)] * spin[(x - 1)*n + ( j + 1)] - svaz[(2 * x - 1)*(n+1) + ( j)] * spin[(x - 2)*n + ( j)] * spin[(x - 1)*n + ( j)] - svaz[(2 * x)*(n+1) + ( j)] * spin[(x - 1)*n + ( j)] * spin[(0)*n + ( j)] - Polevar;
    for (int i = 1; i < x - 1; i++)
        E[i*n + (x - 1)] = -svaz[(i)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i)*n + ( x - 2)] - svaz[(i)*(n+1) + ( x)] * spin[(i)*n + ( x - 1)] * spin[(i)*n + ( 0)] - svaz[(i + x)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i - 1)*n + ( x - 1)] - svaz[(i + x + 1)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i + 1)*n + ( x - 1)] - Polevar;
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E1[i*n +j] = E[i*n +j];
        }
    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(short *E,float *E1)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E1[i*n +j] = E[i*n +j];
        }
    }
}
///////////////////////////////////////////////////////Копирует массив с энергией Е в Е1
void CopyE(float *E,short *E1)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            E1[i*n +j] = E[i*n +j];
        }
    }
}
////////////////////////////////////////////////////
void CopyE(short *E,short *E1)
{
    for (int i = 0; i < n; i++)
    {
        for (int j=0;j<n;j++)
        {
            E1[i*n +j] = E[i*n +j];
        }

    }
}
/////////////////////////////////////////////////////////Подсчитывает энергию одного спина
float EnergyIJ(int i,int j,short spin[][n],short svaz[][n+1])
{
    int x = n;
    float Ener = 0;
    if (i == 0 && j == 0)
        Ener = -svaz[0][ 0] * spin[0][ 0] * spin[0][ x - 1] - svaz[0][ 1] * spin[0][ 0] * spin[0][ 1] - svaz[x][ 0] * spin[0][ 0] * spin[x - 1][ 0] - svaz[x + 1][ 0] * spin[0][ 0] * spin[1][ 0] - Polevar;
    if (i == 0 && j == x - 1)
        Ener = -svaz[0][ x - 1] * spin[0][ x - 1] * spin[0][ x - 2] - svaz[0][ x] * spin[0][ x - 1] * spin[0][ 0] - svaz[x][ x - 1] * spin[0][ x - 1] * spin[x - 1][ x - 1] - svaz[x + 1][ x - 1] * spin[0][ x - 1] * spin[1][ x - 1] - Polevar;
    if (i == x - 1 && j == 0)
        Ener = -svaz[x - 1][ 0] * spin[x - 1][ x - 1] * spin[x - 1][ 0] - svaz[x - 1][ 1] * spin[x - 1][ 0] * spin[x - 1][ 1] - svaz[2 * x - 1][ 0] * spin[x - 2][ 0] * spin[x - 1][ 0] - svaz[2 * x][ 0] * spin[x - 1][ 0] * spin[0][ 0] - Polevar;
    if (i == x - 1 && j == x - 1)
        Ener = -svaz[x - 1][ x - 1] * spin[x - 1][ x - 1] * spin[x - 1][ x - 2] - svaz[x - 1][ x] * spin[x - 1][ x - 1] * spin[x - 1][ 0] - svaz[2 * x - 1][ x - 1] * spin[x - 2][ x - 1] * spin[x - 1][ x - 1] - svaz[2 * x][ x - 1] * spin[x - 1][ x - 1] * spin[0][ x - 1] - Polevar;
    if (i != 0 && j != 0 && i != x - 1 && j != x - 1)
        Ener = -svaz[i][ j] * spin[i][ j - 1] * spin[i][ j] - svaz[i][ j + 1] * spin[i][ j] * spin[i][ j + 1] - svaz[i + x][ j] * spin[i - 1][ j] * spin[i][ j] - svaz[i + x + 1][ j] * spin[i][ j] * spin[i + 1][ j] - Polevar;
    if (i == 0 && j != 0 && j != x - 1)
        Ener = -svaz[0][ j] * spin[0][ j] * spin[0][ j - 1] - svaz[0][ j + 1] * spin[0][ j] * spin[0][ j + 1] - svaz[x][ j] * spin[0][ j] * spin[x - 1][ j] - svaz[x + 1][ j] * spin[0][ j] * spin[1][ j] - Polevar;
    if (j == 0 && i != 0 && i != x - 1)
        Ener = -svaz[i][ 0] * spin[i][ 0] * spin[i][ x - 1] - svaz[i][ 1] * spin[i][ 0] * spin[i][ 1] - svaz[i + x][ 0] * spin[i][ 0] * spin[i - 1][ 0] - svaz[i + x + 1][ 0] * spin[i][ 0] * spin[i + 1][ 0] - Polevar;
    if (i == x - 1 && j != 0 && j != x - 1)
        Ener = -svaz[x - 1][ j] * spin[x - 1][ j - 1] * spin[x - 1][ j] - svaz[x - 1][ j + 1] * spin[x - 1][ j] * spin[x - 1][ j + 1] - svaz[2 * x - 1][ j] * spin[x - 2][ j] * spin[x - 1][ j] - svaz[2 * x][ j] * spin[x - 1][ j] * spin[0][ j] - Polevar;
    if (j == x - 1 && i != 0 && i != x - 1)
        Ener = -svaz[i][ x - 1] * spin[i][ x - 1] * spin[i][ x - 2] - svaz[i][ x] * spin[i][ x - 1] * spin[i][ 0] - svaz[i + x][ x - 1] * spin[i][ x - 1] * spin[i - 1][ x - 1] - svaz[i + x + 1][ x - 1] * spin[i][ x - 1] * spin[i + 1][ x - 1] - Polevar;
    return Ener;
}
/////////////////////////////////////////////////////////Подсчитывает энергию одного спина
float EnergyIJ(int i,int j,short *spin,short *svaz)
{
    int x = n;
    float Ener = 0;
    if (i == 0 && j == 0)
        Ener = -svaz[(0)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(0)*n + ( x - 1)] - svaz[(0)*(n+1) + ( 1)] * spin[(0)*n + ( 0)] * spin[(0)*n + ( 1)] - svaz[(x)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(x - 1)*n + ( 0)] - svaz[(x + 1)*(n+1) + ( 0)] * spin[(0)*n + ( 0)] * spin[(1)*n + ( 0)] - Polevar;
    if (i == 0 && j == x - 1)
        Ener = -svaz[(0)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(0)*n + ( x - 2)] - svaz[(0)*(n+1) + ( x)] * spin[(0)*n + ( x - 1)] * spin[(0)*n + ( 0)] - svaz[(x)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] - svaz[(x + 1)*(n+1) + ( x - 1)] * spin[(0)*n + ( x - 1)] * spin[(1)*n + ( x - 1)] - Polevar;
    if (i == x - 1 && j == 0)
        Ener = -svaz[(x - 1)*(n+1) + ( 0)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( 0)] - svaz[(x - 1)*(n+1) + ( 1)] * spin[(x - 1)*n + ( 0)] * spin[(x - 1)*n + ( 1)] - svaz[(2 * x - 1)*(n+1) + ( 0)] * spin[(x - 2)*n + ( 0)] * spin[(x - 1)*n + ( 0)] - svaz[(2 * x)*(n+1) + ( 0)] * spin[(x - 1)*n + ( 0)] * spin[(0)*n + ( 0)] - Polevar;
    if (i == x - 1 && j == x - 1)
        Ener = -svaz[(x - 1)*(n+1) + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 2)] - svaz[(x - 1)*(n+1) + ( x)] * spin[(x - 1)*n + ( x - 1)] * spin[(x - 1)*n + ( 0)] - svaz[(2 * x - 1)*(n+1) + ( x - 1)] * spin[(x - 2)*n + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] - svaz[(2 * x)*(n+1) + ( x - 1)] * spin[(x - 1)*n + ( x - 1)] * spin[(0)*n + ( x - 1)] - Polevar;
    if (i != 0 && j != 0 && i != x - 1 && j != x - 1)
        Ener = -svaz[(i)*(n+1) + ( j)] * spin[(i)*n + ( j - 1)] * spin[(i)*n + ( j)] - svaz[(i)*(n+1) + ( j + 1)] * spin[(i)*n + ( j)] * spin[(i)*n + ( j + 1)] - svaz[(i + x)*(n+1) + ( j)] * spin[(i - 1)*n + ( j)] * spin[(i)*n + ( j)] - svaz[(i + x + 1)*(n+1) + ( j)] * spin[(i)*n + ( j)] * spin[(i + 1)*n + ( j)] - Polevar;
    if (i == 0 && j != 0 && j != x - 1)
        Ener = -svaz[(0)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(0)*n + ( j - 1)] - svaz[(0)*(n+1) + ( j + 1)] * spin[(0)*n + ( j)] * spin[(0)*n + ( j + 1)] - svaz[(x)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(x - 1)*n + ( j)] - svaz[(x + 1)*(n+1) + ( j)] * spin[(0)*n + ( j)] * spin[(1)*n + ( j)] - Polevar;
    if (j == 0 && i != 0 && i != x - 1)
        Ener = -svaz[(i)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i)*n + ( x - 1)] - svaz[(i)*(n+1) + ( 1)] * spin[(i)*n + ( 0)] * spin[(i)*n + ( 1)] - svaz[(i + x)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i - 1)*n + ( 0)] - svaz[(i + x + 1)*(n+1) + ( 0)] * spin[(i)*n + ( 0)] * spin[(i + 1)*n + ( 0)] - Polevar;
    if (i == x - 1 && j != 0 && j != x - 1)
        Ener = -svaz[(x - 1)*(n+1) + ( j)] * spin[(x - 1)*n + ( j - 1)] * spin[(x - 1)*n + ( j)] - svaz[(x - 1)*(n+1) + ( j + 1)] * spin[(x - 1)*n + ( j)] * spin[(x - 1)*n + ( j + 1)] - svaz[(2 * x - 1)*(n+1) + ( j)] * spin[(x - 2)*n + ( j)] * spin[(x - 1)*n + ( j)] - svaz[(2 * x)*(n+1) + ( j)] * spin[(x - 1)*n + ( j)] * spin[(0)*n + ( j)] - Polevar;
    if (j == x - 1 && i != 0 && i != x - 1)
        Ener = -svaz[(i)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i)*n + ( x - 2)] - svaz[(i)*(n+1) + ( x)] * spin[(i)*n + ( x - 1)] * spin[(i)*n + ( 0)] - svaz[(i + x)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i - 1)*n + ( x - 1)] - svaz[(i + x + 1)*(n+1) + ( x - 1)] * spin[(i)*n + ( x - 1)] * spin[(i + 1)*n + ( x - 1)] - Polevar;
    return Ener;
}
///////////////////////////////////////////////////////////
void MCProhod(double temp,short *spin,short *svaz,float *E)
{

    float En1,En2;
    int step=0;
    double slch;
    double veroyatnost=0;
    while(step<=prohod_MC2)
    {
        int slspin=rand()%(n*n);
        int i = slspin / n;
        int j = slspin - i * n;
        int var = spin[i*n +j];
        En1 = EnergyIJ(i, j, spin, svaz);
        spin[i*n + j] *= -1;
        En2 = EnergyIJ(i, j, spin, svaz);
        slch=((double)rand())/RAND_MAX;
        if (En2 <= En1)
            veroyatnost = 1;
        else
            veroyatnost = exp(-(En2 - En1) / temp);
        if (slch < veroyatnost){
            E[i*n +j] = En2;
            //sva
            if(i==0){
                E[(n-1)*n+j]=EnergyIJ((n-1), j, spin, svaz);
            }
            else{
                E[slspin-n]=EnergyIJ(i-1, j, spin, svaz);
            }
            //sn
            if(i==n-1){
                E[j]=EnergyIJ(0, j, spin, svaz);
            }
            else{
                E[slspin+n]=EnergyIJ(i+1, j, spin, svaz);
            }
            //sl
            if(j==0){
                E[i*n + (n-1)]=EnergyIJ(i, (n-1), spin, svaz);
            }
            else{
                E[slspin - 1]=EnergyIJ(i, j-1, spin, svaz);
            }
            //sp
            if(j==n-1){
                E[i*n]=EnergyIJ(i,0, spin, svaz);
            }
            else{
                E[slspin + 1]=EnergyIJ(i, j+1, spin, svaz);
            }
        }
        else
            spin[i*n +j] *= -1;
        step++;
    }
    //Energy(spin,svaz,E);
}
///////////////////////////////////////////////////////////Подсчет намагниченности
double MagnMet(short *spin)
{
    double Magnvar=0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Magnvar += spin[i*n +j];
        }
    }
    return abs(Magnvar/(n*n));
}
////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int Claster(float *E1)
{
    int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (E1[i*n +j] == -4 || E1[i*n +j] == -2)
            {
                tmp = 1;
                int top = i * n + j;
                t = MaxClass(top, E1, tmp,Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }
        }
    }
    delete []Ochered;
    return Max;
}
//////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int ClasterFerr(float *E1)
{
    int tmp = 0;
    int t;
    int Max = 0;
    int *Ochered;
    Ochered = new int[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            if (E1[i*n +j] == -4 )
            {
                tmp = 1;
                int top = i * n + j;
                t = MaxClassFerr(top, E1, tmp, Ochered);
                if (t >= Max)
                    Max = t;
                t = 0;
            }
        }
    }
    delete []Ochered;
    return Max;
} 
//////////////////////////////////////////////////////Для подсчета макс.кластера для спинового стекла
int MaxClass(int per,float *E1,int tmp,int *Ochered)
{
    int sl, sp, sv, sn, top, t1, t2;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    t1 = per % n;
    t2 = (per - t1) / n;
    E1[t2*n +t1] = 0;
    while (w < r)
    {
        top = Ochered[w];
        ////
        sv = sosed_sv(top);
        t1 = sv % n;
        t2 = (sv - t1) / n;
        if (E1[t2*n +t1] == -4 || E1[t2*n +t1] == -2)
        {
            Ochered[r] = sv;
            E1[t2*n +t1] = 0;
            r++;
        }
        ////
        sn = sosed_sn(top);
        t1 = sn % n;
        t2 = (sn - t1) / n;
        if (E1[t2*n +t1] == -4 || E1[t2*n +t1] == -2)
        {
            Ochered[r] = sn;
            E1[t2*n +t1] = 0;
            r++;
        }
        ////
        sl = sosed_sl(top);
        t1 = sl % n;
        t2 = (sl - t1) / n;
        if (E1[t2*n +t1] == -4 || E1[t2*n +t1] == -2)
        {
            Ochered[r] = sl;
            E1[t2*n +t1] = 0;
            r++;
        }
        ////
        sp = sosed_sp(top);
        t1 = sp % n;
        t2 = (sp - t1) / n;
        if (E1[t2*n +t1] == -4 || E1[t2*n +t1] == -2)
        {
            Ochered[r] = sp;
            E1[t2*n +t1] = 0;
            r++;
        }
        ////
        w++;

    }
    return r;


}
/////////////////////////////////////////////////////Для подсчета макс.кластера для ферромагеника
int MaxClassFerr(int per,float *E1,int tmp, int *Ochered)
{
    int sl, sp, sv, sn, top, t1, t2;
    int w = 0;
    int r = 1;
    Ochered[0] = per;
    t1 = per % n;
    t2 = (per - t1) / n;
    E1[t2*n + t1] = 0;
    while (w < r)
    {
        top = Ochered[w];
        ////
        sv = sosed_sv(top);
        t1 = sv % n;
        t2 = (sv - t1) / n;
        if (E1[t2*n + t1] == -4 )
        {
            Ochered[r] = sv;
            E1[t2*n + t1] = 0;
            r++;
        }
        ////
        sn = sosed_sn(top);
        t1 = sn % n;
        t2 = (sn - t1) / n;
        if (E1[t2*n + t1] == -4 )
        {
            Ochered[r] = sn;
            E1[t2*n + t1] = 0;
            r++;
        }
        ////
        sl = sosed_sl(top);
        t1 = sl % n;
        t2 = (sl - t1) / n;
        if (E1[t2*n + t1] == -4 )
        {
            Ochered[r] = sl;
            E1[t2*n + t1] = 0;
            r++;
        }
        ////
        sp = sosed_sp(top);
        t1 = sp % n;
        t2 = (sp - t1) / n;
        if (E1[t2*n + t1] == -4 )
        {
            Ochered[r] = sp;
            E1[t2*n + t1] = 0;
            r++;
        }
        ////
        w++;

    }
    return r;
}
//////////////////////////////////////////////////// Служебные функции для посчета мах.кластера
int sosed_sv( int top)
{
    int j = top % n;
    int i = (top - j) / n;
    int sv = 0;
    if (i == 0)
        sv = n * (n - 1) + j;
    else
        sv = (i - 1) * n + j;
    return sv;
}
int sosed_sn( int top)
{
    int j = top % n;
    int i = (top - j) / n;
    int sn = 0;
    if (i == n-1)
        sn = j;
    else
        sn = (i + 1) * n + j;
    return sn;
}
int sosed_sl( int top)
{
    int j = top % n;
    int i = (top - j) / n;
    int sl = 0;
    if (j == 0)
        sl = i * n + (n - 1);
    else
        sl = i * n + (j - 1);
    return sl;
}
int sosed_sp( int top)
{
    int j = top % n;
    int i = (top - j) / n;
    int sp = 0;
    if (j == n - 1)
        sp = i * n;
    else
        sp = i * n + (j + 1);
    return sp;
}
//////////////////////////////////////////////////////////Суммирует энергию системы
float SumEnergy( float *E)
{
    float Sum = 0;
    for (int i = 0; i < n*n; i++)
    {
        Sum += E[i];
    }
    return Sum/2;
}
/////////////////////////////////////////////////////////// Получает поле

int main(int argc, char **argv)
{
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    /////////////////// время///////////////
    srand((unsigned)time(NULL)+rank);
    //srand(1);
    ////////////////////////////////////////


    int CO = Prohod_MC_sampling/CountOut;
    /////////////////////

    ///////// Создание массива spin[n][n]////////////
    short *spin=new short[n*n];
    /////////////////////////////////////   
    //Создание массива svaz][(n * 2 + 1)][ n + 1]////////////
    short *svaz = new short[(n*2+1)*(n+1)];
    /////////////////////////////////////////

    ////// Создание массива E[n][ n] ////////////

    float *E=new float[n*n];
    /////////////////////////////////////

    ////// Создание массива E1[n][ n] ////////////
    float *E1=new float[n*n];
    ////////////////////////////////////

    ////////// Задание связей ферромагнетика///////////
    // for (int i=0;i<(n*2+1)*(n+1);i++)
    // {
    //     svaz[i]=1;

    // }
    /////////////массивы для вывода//////////////////////////

    //if(current_value)//
//    double *recvE = new double[CountOut*size];
//    double *recvM = new double[CountOut*size];
//    double *recvPP = new double[CountOut*size];
//    double *recvPPF = new double[CountOut*size];
    /////////////////////


    double *recvAE = new double[size];
    double *recvAE2 = new double[size];
    double *recvAE4 = new double[size];
    double *recvAM = new double[size];
    double *recvAM2 = new double[size];
    double *recvAM4 = new double[size];
    double *recvAPP = new double[size];
    double *recvAPPF = new double[size];

    double *recvHeatCapacity = new double[size];
    double *recvVospr = new double[size];

    double *recvBCenergy = new double[size];
    double *recvBCmagn = new double[size];
    double *recvBCPPF = new double[size];
    double *recvBCPP = new double[size];

//////////////////////////
//    //if(current_value)
//    vector<double> CValue;          // E
//    vector<double> CValuePP;        // PP
//    vector<double> CValuePPF;       // PPF
//    vector<double> CValueMagn;      // M

//    ofstream outE("energy.dat",ios::out);
//    ofstream outPP("PP.dat",ios::out);
//    ofstream outM("Magn.dat",ios::out);
//    ofstream outPPF("PPF.dat",ios::out);
/////////////////////////


    ofstream outAPP("APP.dat",ios::out);            //средний параметр порядка
    ofstream outAPPF("APPF.dat",ios::out);          //средний параметр порядка2
    ofstream outAE("energyAVG.dat",ios::out);       //средняя энергия
    ofstream outAE2("energyAVG2.dat",ios::out);     //средняя энергия квадрат
    ofstream outAE4("energyAVG4.dat",ios::out);     //средняя энергия 4 степень
    ofstream outAM("AMagn.dat",ios::out);           //Средняя намагниченность
    ofstream outHeatCapacity("C.dat",ios::out);     //Теплоемкость
    ofstream outAM2("AMagn2.dat",ios::out);			
    ofstream outAM4("AMagn4.dat",ios::out);
    ofstream outVospr("X.dat",ios::out);            //Восприимчивость


    // биндеры
    ofstream outBCenergy("BC_energy.dat",ios::out); //биндер по энергии
    ofstream outBCmagn("BC_magn.dat",ios::out);     //биндер по намагниченности
    ofstream outBCPP("BC_PP.dat",ios::out);         //биндер по параметру порядка
    ofstream outBCPPF("BC_PPF.dat",ios::out);       //биндер по параметру порядка2
    ////////////////////////////
    double error_svaz;

    // создание связей с проверкой
    //  error_svaz=Slsvaz(svaz,num_podhod); 
    // cout<<"Kolichestvo oshibok: "<<error_svaz<<"  ot "<<rank<<'\n';
    Randsvaz(svaz); // случайное создание связей
    ////////////////////////////
    SpinFerr(spin);   // упорядоченное состояние
    //Slspinn(spin);      // случайное состояние
    Energy(spin,svaz,E);
    CopyE(E,E1);
    // Vivoderrors(n,svaz);
    //if(rank==0)
        //VivodSSS(n,spin,svaz);
    //////MC///////////////////////////////


    double SumenergyVar,MagnVar,MagnVar2,MagnVar4;
    double maxSF,PPF,PPF2,PPF4,maxS,PP,PP2,PP4;

    double temPP;
    double temPPF;
    double temMagn;

    double SumenergyVar2, SumenergyVar4;
    double Aenergy=0;
    double Aenergy2=0;
    double Aenegry4=0;
    double Amagn2=0;
    double Amagn4=0;
    double APPF2=0;
    double APPF4=0;
    double APP2=0;
    double APP4=0;
    double X_o=0;

    //++++++++++++++++++++++++++++++++++++++
    double e1=0;

    // для вывода БК for PP2

    double RBCenergy=0;
    double RBCmagn=0;
    double RBCPPF=0;
    double RBCPP=0;


    // // для вывода C
    double C_o=0;


    //++++++++++++++++++++++++++++++++++++++
    //для репличного обмена

    unsigned long long index=0;
    double exchange_t=0;
    double exchange_e=0;
    int yes_no_exchange=0;

    double probability_of_exchange=0;
    MPI_Status status;
    short *exchange_state=new short[n*n];

    double temp=(rank+1)*(maxtemp-mintemp)/size+mintemp;
    unsigned long long Prohod=0;
    int schetchik=0;
    cout<<"hello from processor"<<rank<<endl;
    double t1, t2;
    double prev_time=0; 
	t1 = MPI_Wtime();
    ////////////////////////Разогрев
    while(Prohod<Prohod_MC_equilibration)
    {
        schetchik++;
        //Energy(spin,svaz,E);
        e1=SumEnergy(E);
        
        if(replica_on)
        {
            if(Prohod%between_exchange_equilibration==0)
            {
                for(int i=size-1;i>0;--i)
                {
                    if(rank==i)
                    {
                        MPI_Send(&temp,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Send(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&yes_no_exchange,1, MPI_INT, i-1, 0, MPI_COMM_WORLD,&status);
                        if(!yes_no_exchange);

                        else
                        {
                            MPI_Send(spin,n*n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD);
                            MPI_Recv(exchange_state,n*n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD,&status);
                            CopyE(exchange_state,spin);
                        }
                    }
                    if(rank==i-1)
                    {
                        MPI_Recv(&exchange_t,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,&status);
                        MPI_Recv(&exchange_e,1, MPI_DOUBLE, i, 0,MPI_COMM_WORLD,&status);
                        probability_of_exchange=exp(((1/temp)-(1/exchange_t))*(e1-exchange_e));
                        if (probability_of_exchange<(double)(rand())/RAND_MAX)
                        {
                            yes_no_exchange=0;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD);

                        }
                        else
                        {
                            yes_no_exchange=1;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD); //порядок?
                            MPI_Recv(exchange_state,n*n, MPI_SHORT, i, 0, MPI_COMM_WORLD,&status);
                            MPI_Send(spin,n*n, MPI_SHORT, i, 0, MPI_COMM_WORLD);
                            MPI_Send(&e1,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD); //порядок?
                            e1=exchange_e;
                            CopyE(exchange_state, spin);
                        }
                    }

                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        Energy(spin,svaz,E);
        MCProhod(temp,spin,svaz,E);
        if(rank==0)
        {
            if(Prohod%(Prohod_MC_equilibration/100)==0) //было if(Prohod%(Prohod_MC_sampling/100)==0)
            {
            	t2=MPI_Wtime();
                cout<<"Equilibration status: "<<(Prohod)/(Prohod_MC_equilibration/100)+1<<" time per status: "<<t2-t1-prev_time<<endl; // было cout<<"Status: "<<Prohod/(Prohod_MC_sampling/100)<<endl;
                prev_time=t2-t1;
            }    
        }
        Prohod++;
    }
    Prohod=0;
    prev_time=0;
    while(Prohod<Prohod_MC_sampling)
    {
        //Energy(spin,svaz,E);
        e1=SumEnergy(E);
        if(replica_on)
        {
            if(Prohod%between_exchange==0)
            {
                for(int i=size-1;i>0;--i)
                {
                    if(rank==i)
                    {
                        MPI_Send(&temp,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Send(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&yes_no_exchange,1, MPI_INT, i-1, 0, MPI_COMM_WORLD,&status);
                        if(!yes_no_exchange);

                        else
                        {
                            MPI_Send(spin,n*n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD);
                            MPI_Recv(exchange_state,n*n, MPI_SHORT, i-1, 0, MPI_COMM_WORLD,&status);
                            MPI_Recv(&e1,1, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD,&status);
                            CopyE(exchange_state,spin);
                        }
                    }
                    if(rank==i-1)
                    {
                        MPI_Recv(&exchange_t,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,&status);
                        MPI_Recv(&exchange_e,1, MPI_DOUBLE, i, 0,MPI_COMM_WORLD,&status);
                        probability_of_exchange=exp(((1/temp)-(1/exchange_t))*(e1-exchange_e));
                        double ttt1=(double)(rand())/RAND_MAX;
                        if (probability_of_exchange<ttt1)
                        {
                            yes_no_exchange=0;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD);

                        }
                        else
                        {
                            yes_no_exchange=1;
                            MPI_Send(&yes_no_exchange,1, MPI_INT, i, 0, MPI_COMM_WORLD); //порядок?
                            MPI_Recv(exchange_state,n*n, MPI_SHORT, i, 0, MPI_COMM_WORLD,&status);
                            MPI_Send(spin,n*n, MPI_SHORT, i, 0, MPI_COMM_WORLD);
                            MPI_Send(&e1,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD); //порядок?
                            CopyE(exchange_state, spin);
                            e1=exchange_e;

                        }
                    }

                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }

        Energy(spin,svaz,E);
        MCProhod(temp,spin,svaz,E);

        if(Prohod%CO==0)
        {
            index++;
            if (Polevar==0)
            {
                CopyE(E,E1);                    // при расчете кластера, коцаеца массив энергий, поэтому мы его копируем и работаем с копией
                // maxSF = ClasterFerr(E1);
                // PPF = maxSF / (n*n);                  //
                // PPF2 = PPF*PPF;
                // PPF4 = PPF*PPF*PPF*PPF;

                // temPPF=(PPF+(index-1)*temPPF)/(index);      // не забыть помножить !!!*N куммулятив ферромагнитного кластера
                // APPF2=(PPF2+(index-1)*APPF2)/(index);
                // APPF4=(PPF4+(index-1)*APPF4)/(index);

                // cout<<"OK1"<<endl;

                CopyE(E,E1);                        // при расчете кластера, коцаеца массив энергий, поэтому мы его копируем и работаем с копией
                maxS = Claster(E1);
                PP= maxS / (n*n);
                PP2 = PP*PP;
                PP4 = PP*PP*PP*PP;
                temPP=(PP+(index-1)*temPP)/(index);      // !!!*N куммулятив смешанного кластера
                APP2=(PP2+(index-1)*APP2)/(index);
                APP4=(PP4+(index-1)*APP4)/(index);
            }
            // cout<<"OK2"<<endl;

            SumenergyVar = SumEnergy(E);     /////////////
            SumenergyVar2=pow(SumenergyVar,2);
            SumenergyVar4=pow(SumenergyVar,4);
            Aenergy=(SumenergyVar+(index-1)*Aenergy)/(index);   //было tem213=Avalue*(CO)/(Prohod);
            Aenergy2=(SumenergyVar2+(index-1)*Aenergy2)/(index);
            Aenegry4=(SumenergyVar4+(index-1)*Aenegry4)/(index);

            // cout<<"OK3"<<endl;

            MagnVar=MagnMet(spin);			/////////////
            MagnVar2=MagnVar*MagnVar;
            MagnVar4=MagnVar*MagnVar*MagnVar*MagnVar;
            temMagn=(MagnVar+(index-1)*temMagn)/(index);   // !!!*N куммулятив смешанного кластера
            Amagn2=(MagnVar2+(index-1)*Amagn2)/(index);
            Amagn4=(MagnVar4+(index-1)*Amagn4)/(index);

            // cout<<"OK4"<<endl;

//            if(current_value){
//                CValue.insert(CValue.end(),SumenergyVar);       // E
//                CValueMagn.insert(CValueMagn.end(),MagnVar);    // M
//                CValuePP.insert(CValuePP.end(),PP);             // PP
//                CValuePPF.insert(CValuePPF.end(),PPF);          //PPF
//            }


            if(rank==0)
            {
                if(Prohod%(Prohod_MC_sampling/100)==0) //было if(Prohod%(Prohod_MC_sampling/100)==0)
                {
                	t2=MPI_Wtime();
                    cout<<"Status: "<<(Prohod)/(Prohod_MC_sampling/100)+1<<" time per status: "<<t2-t1-prev_time<<endl; // было cout<<"Status: "<<Prohod/(Prohod_MC_sampling/100)<<endl;
                    prev_time=t2-t1;
                }
            }
        }
        Prohod++;
    }

    //cout<<"index:"<<index<<"    Coutout:"<<CountOut<<endl;


    C_o=((Aenergy2)-(Aenergy*Aenergy))/(temp*temp*n*n); //  теплоёмкость //  теплоёмкость
    X_o=((Amagn2)-(temMagn*temMagn))/(temp*n*n);

    //кумулянт Биндера для энергии
    RBCenergy = 1 -(Aenegry4/(3*pow(Aenergy2,2)));  // результирующий БК энергии
    RBCmagn = 1-(Amagn4/(3*pow(Amagn2,2)));  // результирующий БК магнитный
    RBCPPF = 1-(APPF4/(3*pow(APPF2,2)));  // результирующий БК параметра порядка 2
    RBCPP = 1-(APP4/(3*pow(APP2,2)));  // результирующий БК параметра порядка

    MPI_Gather(&X_o, 1, MPI_DOUBLE, recvVospr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&C_o, 1, MPI_DOUBLE, recvHeatCapacity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Теплоемкость
    MPI_Gather(&RBCenergy, 1, MPI_DOUBLE, recvBCenergy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК энергию
    MPI_Gather(&RBCmagn, 1, MPI_DOUBLE, recvBCmagn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК магн
    MPI_Gather(&RBCPPF, 1, MPI_DOUBLE, recvBCPPF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PPF
    MPI_Gather(&RBCPP, 1, MPI_DOUBLE, recvBCPP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // БК PP
    MPI_Gather(&Aenergy, 1, MPI_DOUBLE, recvAE, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenergy2, 1, MPI_DOUBLE, recvAE2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Aenegry4, 1, MPI_DOUBLE, recvAE4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temMagn, 1, MPI_DOUBLE, recvAM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Gather(&Amagn2, 1, MPI_DOUBLE, recvAM2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&Amagn4, 1, MPI_DOUBLE, recvAM4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temPP, 1, MPI_DOUBLE, recvAPP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&temPPF, 1, MPI_DOUBLE, recvAPPF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//    if(current_value){
//        MPI_Gather(&CValue[0], CountOut, MPI_DOUBLE, recvE, CountOut, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        MPI_Gather(&CValueMagn[0], CountOut, MPI_DOUBLE, recvM, CountOut, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        MPI_Gather(&CValuePP[0], CountOut, MPI_DOUBLE, recvPP, CountOut, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        MPI_Gather(&CValuePPF[0], CountOut, MPI_DOUBLE, recvPPF, CountOut, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    }

    if(rank==0) //countOut - количество выводов!!! CO- какое значение из общего количества
    {
        for(int i=0; i<size; ++i)
        {

//            if(current_value){
//                for(int j=0; j<CountOut; ++j)
//                {
//                    outE<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvE[i*CountOut+j]<<endl;
//                    outPP<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvPP[i*CountOut+j]<<endl;
//                    outPPF<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvPPF[i*CountOut+j]<<endl;
//                    outM<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvM[i]<<endl;
//                }
//            }
         //   outAPPF<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPPF[i]<<endl;
            outAPP<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAPP[i]<<endl;
            outAM<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM[i]<<endl;
             outAM2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM2[i]<<endl;
              outAM4<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAM4[i]<<endl;
            outAE<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE[i]<<endl;
            outAE2<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE2[i]<<endl;
            outAE4<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvAE4[i]<<endl;

            outHeatCapacity<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvHeatCapacity[i]<<endl;                 // Теплоемкость

            outVospr<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvVospr[i]<<endl;  //воспириимчивость
            outBCenergy<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCenergy[i]<<endl; // БК energy
            outBCmagn<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCmagn[i]<<endl; // БК magn
           // outBCPPF<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPPF[i]<<endl; // БК PPF
            outBCPP<<(i+1)*((maxtemp-mintemp)/size)+mintemp<<'\t'<<recvBCPP[i]<<endl; // БК PP
        }
    }
    t2=MPI_Wtime();
    cout<<"total time: "<<t2-t1<<endl;

//    //if(current_value)
//    delete []recvE;
//    delete []recvM;
//    delete []recvPP;
//    delete []recvPPF;

    delete []svaz;
    delete []E;
    delete []E1;
    delete []spin;
    delete []exchange_state;
    delete []recvAPPF;
    delete []recvAPP;
    delete []recvAM;
    delete []recvAM2;
    delete []recvAM4;
    delete []recvAE;
    delete []recvAE2;
    delete []recvAE4;
    delete []recvBCenergy;
    delete []recvBCmagn;
    delete []recvBCPPF;
    delete []recvBCPP;
    delete []recvHeatCapacity;
    delete []recvVospr;
    MPI_Finalize();
    cout<<"FINISH"<<endl;


    return 0;
}


