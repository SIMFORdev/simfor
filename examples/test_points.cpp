#include <stdlib.h>
#include <iostream>
#include <GL/glut.h>
#include "../include/simfor/plotter.h"


using namespace std;


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
    Show_Legend(false);              //Баг с масштабированием

    glutInit(&argc, argv); // Инициализация FreeGLUT             
	
    // Легенда
    std::string legend_text[1];
	legend_text[0] = "points";
    Set_Legend_Text(legend_text);    

    // Размер окна
    initialization_general(800, 800);

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

	// Первый индекс - ось (0 - x, 1 - y, 2 - z)
	// Второй индекс - номер точки
	// Третий индекс - номер набора
	points[0][0][0] = -3.2;
	points[0][1][0] = -2.3;
	points[0][2][0] = -0.6;
	points[0][3][0] = 1.7;
	points[0][4][0] = 3.3;

	points[1][0][0] = 2;
	points[1][1][0] = 3.5;
	points[1][2][0] = 3;
	points[1][3][0] = 1.7;
	points[1][4][0] = 2.2;

    Set_Graphic_Point_Size(10);
    color_select("Blue");

    // Сейчас нельзя переключаться между набором точек и наборо функций
    // TODO Добавить переключение между набором точек и функций
	initialization_p_many(points, n, 2, 1); 

    OpenGL();

    return 0;
}