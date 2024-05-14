#include <stdlib.h>
#include <iostream>
#include <GL/glut.h>
#include "simfor/plotter.h"


using namespace std;


double sine(double x)
{
    return sin(x);
}

double cosine(double x)
{
    return cos(x);
}

double x2(double x)
{
	return x*x;
}


int main(int argc, char** argv)
{
    // Начальные настройки
    //Set_Axis_X_Text("Ось x");
    //Set_Axis_Y_Text("Ось y");
    //Set_Axis_Z_Text("Ось z");

    //Set_Axis_Text_Size(24.);
    //Set_Text_Line_Size(100.);
    //Set_Marks_Text_Size(24.);
    //Show_Z_Marks(true);

    //Set_Axis_Point_Size(24.);       //Работает
    //Set_Axis_Line_Size(12.);        //Работает
    //Set_Marks_Line_Size(10.);       //Работает

    // График
    //surface_mode(true);
    //Set_Step_X(1);                  //Работает           
    //Set_Step_Y(1);                  //Работает
    //Set_Graphic_Line_Size(6.);      //Работает
    //Set_Graphic_Point_Size(24.);    //Работает

    //Сетка
    //Set_Net(true);                  //Работает                   
    //Set_Net_Line_Size(10);          //Работает
    //tenth_mode(true);               //Работает  

    //Set_Legend_Start_Pos_X(-0.01);  
    Show_Legend(false);              // TODO Баг с масштабированием

    glutInit(&argc, argv); // Инициализация FreeGLUT             
	
    // Легенда
    std::string legend_text[3];
	legend_text[0] = "sin x";
	legend_text[1] = "cos x";
    legend_text[2] = "x^2";
    Set_Legend_Text(legend_text);    

    // Размер окна
    initialization_general(800, 800);

    // Набор функций
    // Сейчас интервал отображения определяется
    // по последней инициализированной функции

    // TODO Вычислять интервал отображения исходя из
    // всех интервалов
	double (*functions_arr[3])(double) = {sine, cosine, x2};
	initialization_function(-5, 5, functions_arr, 1);
	initialization_function(-5, 5, functions_arr, 2);
    initialization_function(-5, 5, functions_arr, 3);

    //Набор точек
    int n = 5;
	double*** points = new double**[2];
    
	for (int i = 0; i < 2; i++)
	{
        points[i] = new double*[n];
        for (int j = 0; j < n; j++) {
            points[i][j] = new double[2];
        }
    }

    OpenGL();

    return 0;
}