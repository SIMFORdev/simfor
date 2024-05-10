#pragma once
#include <stdlib.h>
#include <iostream>
#include <GL/freeglut.h>
#include <cmath>
#include <string>
#include <algorithm>


// ------------------------------ НАСТРОЙКИ ------------------------------

//ЭКРАН - передаются при инициализации
static int window_width,          //Ширина
           window_height;         //Высота
static double aspect;             //Отношение сторон


// ЦВЕТА по умолчанию
static double color_axis[3];      //Ось
static double color_text[3];      //Текст
static double color_array[3];     //Динамический цвет
static int *color_table;          //Список цветов по умолчанию
static double *surface_colormap;

//Здесь хранятся информация о осях для каждой поверхности
//- для корректной закраски
static double **surface_start_end;

//Флаги освещения
static bool lighting_activate = false,
            lighting_deactivate = false,
			lighting_status = false;

//Тип окрашивания поверхности
static bool color_range = true;


//РАЗМЕРЫ
static double axis_point_size = 1.0;      //Точек оси
static double axis_line_size = 3.0;       //Ось
static double net_line_size = 0.5;        //Сетка

static double graphic_point_size = 5.0;   //Точек графика
static double graphic_line_size = 2.0;    //Линий графика

static double surface_point_size = 5.;    //Поверхность
static double surface_line_size = 5.0;    //Сеть поверхности
static double surface_triangle_volume = 0.5; //Толщина поверхности

static double axis_text_size = 0.002,     //Размер надписей
              marks_text_size = 0.001;    //Размер засечек
static double axis_mark_line_size = 2.0;  //Размер засечек
static double text_size = 1.0;            //Толщина линии текста
static double text_point_size = 0.1;      //Размер точек текста
static double def_text_size_big = 0.0075,
              def_text_size_normal = 0.005,
	          def_text_size_small = 0.0025,
	          def_text_size_very_small = 0.0015;
static double def_point_size_big = 5.,
              def_point_size_normal = 2.5,
	          def_point_size_small = 1.;
static double def_line_thick = 2.0,
              def_line_normal = 1.0,
	          def_line_thin = 0.5;



//ТЕКСТ
static std::string axis_ox_text = "",   //Ось х  
                   axis_oy_text = "",   //Ось y
			       axis_oz_text = "";   //Ось z 
static std::string *object_name;        //Надписи легенды.
static std::string color_name;


//ПАРАМЕТРЫ ИЗОБРАЖЕНИЯ
static bool simplify3d = false;         //Только оси
static bool legend_need = true;        //Нужна ли легенда?
static bool net_need = false;           //Отображать сетку?
static bool mm = false;                 //Миллиметровый режим
static bool antialiasing = false;       //Сглаживание
static double zero_x = 0.66;            //С какой позиции идёт легенда - от - 1 до 1
static bool axis_text_need = true;      //Показ подписей осей
static double diapason_x_start,         //Начало и конец диапазонов отображения
              diapason_x_end,   
              diapason_y_start,
			  diapason_y_end,
			  diapason_z_start,
			  diapason_z_end;
static bool build_z_marks = false;      //Строить засечки для оси Z ?
static double diapason_width,           //Ширина
              diapason_height,          //Высота
	          diapason_depth;           //Глубина


//ВЫЧИСЛЕНИЯ (передаются при инициализации в разных видах)
static int point_number_x,           //Количество точек
           point_number_y,
		   point_number_z;

static int number_of_sets,           //Число наборов объектов
           number_of_intervals = 0;  //Число интервалов
static double step_x = 0.15,          //Шаг по осям
              step_y = 0.15;


//ОТОБРАЖЕНИЕ
static int current_set = 0;          //Текущий отображаемый набор
static bool all_sets = false;        //Изображать все наборы
static double x_x_position = 0,      //Относительное позиционирование
              y_y_position = 0,
			  z_z_position = 0;

//Рисуем - поверхность?
static bool is_surface = false;

//Тип графика линий
//1 - точечный
//2 - линии
//3 - линии с точками
static int line_type = 2;

//Тип поверхности
//1 - точки
//2 - линии и точки
//3 - линии
//4 - поверхность
static int surface_type = 1;


//КАМЕРА
static int window_zoom = 0,
           x_shift = 0,       //Смещение по осям
		   y_shift = 0,
		   z_shift = 0;
static int add_x_canvas = 0,  //Пределы изображения по осям
           add_y_canvas = 0,
		   add_z_canvas = 0;
static double factor = 1,     //Приближение/Удаление от точки
              new_factor = 1;
static double x_degree = 0,   //Вращение
              y_degree = 0,
			  z_degree = 0;


static double *surface_tesseract;
static double(**input_f)(double);

// -------------------------------------------------------------------------------


class Point
{
public:
	bool valid;
	double x, y, z;
	double static_color[3];
	double dynamic_color[3];

	Point()
	{
		valid = true;
		x = 0;
		y = 0;
		z = 0;
		for (int i = 0; i < 3; i++)
		{
			static_color[i] = 0;
			dynamic_color[i] = 0;
		}
	}
	
};
static Point *points_array;


class bound
{
public:
	Point left, right;

	bool eqv(int coord)
	{
		double temp;
		switch (coord)
		{
		case 0:
			temp = left.x - right.x;
			break;
		case 1:
			temp = left.y - right.y;
			break;
		case 2:
			temp = left.z - right.z;
			break;
		default:
			temp = 1;
		}
		temp = abs(temp);
		return temp < 0.00000001;
	}
};
static bound **exclusion_intervals;


bool bounds_check(bound interval, double evalue, int coord);
bool infinity_check(double evalue);
bool nan_check(double evalue);


class Graphic
{
public:
	Point *points;
	double(*fun_y)(double);
	double(*fun_z)(double);
	
	Graphic(double(*new_fun)(double))
	{
		fun_y = new_fun;
		fun_z = new_fun;
		points = new Point[point_number_x + 1];
		Comp();
	}

	Graphic()
	{
		points = new Point[point_number_x + 1];
		fun_y = nullptr;
		fun_z = nullptr;
	}

	//Расчёт точек поверхности
	void Comp()
	{
		for (int i = 0; i < point_number_x; i++)
		{
			bool skip = false;
			points[i].x = 0;
			points[i].y = 0;
			points[i].z = 0;

			double temp = diapason_x_start + step_x*i;
			
			for (int k = 0; k < number_of_intervals; k++)
				if (bounds_check(exclusion_intervals[current_set][k], temp, 0))
					skip = true;
			
			if (infinity_check(temp) || nan_check(temp))
				skip = true;

			if (skip == false)
			{
				points[i].x = temp;
				temp = fun_y(points[i].x);
				for (int k = 0; k < number_of_intervals; k++)
					if (bounds_check(exclusion_intervals[current_set][k], temp, 1))
						skip = true;
				if (infinity_check(temp) || nan_check(temp))
					skip = true;

			}
			if (skip == false)
			{
				points[i].y = temp;
				if (fun_z != nullptr)
				{
					temp = fun_z(points[i].x);
					for (int k = 0; k < number_of_intervals; k++)
						if (bounds_check(exclusion_intervals[current_set][k], temp, 2))
							skip = true;
					if (infinity_check(temp) || nan_check(temp))
						skip = true;
					points[i].z = temp;
				}
			}

			if (skip != false)	points[i].valid = false;
		}
	}
};
static Graphic *function_array;


class Surface
{
public:
	//Поверхность
	double(*fun_surface)(double, double);
	
	Point **point;
	
	double surface_start, surface_end;
	
	Surface(double(*new_f)(double, double))
	{
		fun_surface = new_f;
		point = new Point*[point_number_x];
		for (int i = 0; i < point_number_x; i++)
			point[i] = new Point[point_number_y];

	}

	Surface()
	{
		point = new Point*[point_number_x];
		for (int i = 0; i < point_number_x; i++)
			point[i] = new Point[point_number_y]();
		fun_surface = nullptr;
	}

	void compute_surface()
	{
		surface_tesseract = new double[number_of_sets];

		double min_v = 1000000;
		double max_v = -1000000;

		for (int i = 0; i < point_number_x; i++)
			for (int j = 0; j < point_number_y; j++)
			{
				point[i][j].x = diapason_x_start + step_x*j;
				point[i][j].y = diapason_y_start + step_y*i;
				
				double temp = fun_surface(point[i][j].x, point[i][j].y);
				
				if (infinity_check(temp) || nan_check(temp))
				{
					point[i][j].z = 0;
					point[i][j].valid = true;
				}
				else
				{
					point[i][j].z = -1.*temp;
					if (temp > max_v)
						max_v = temp;
					if (temp < min_v)
						min_v = temp;
				}
			}

		surface_start = min_v;
		surface_end = max_v;

		diapason_z_start = std::min(min_v, diapason_z_start);
		diapason_z_start = floor(diapason_z_start);
		diapason_z_end = std::max(max_v, diapason_z_end);
		diapason_z_end = ceil(diapason_z_end);
	}
};
static Surface *surf_array;


//Вращение, смещение, факторизация
void factorization();

//Получить цвет по номеру
std::string color_map(int numer);
//Получение цвета по названию
void color_select(std::string color_name);
//Фон
void background();
//Общая инициализация, передаются размеры окна
void initialization_general(int i_window_width, int i_window_height);

// Параметры линии
void line_settings(double line_width, double color_r, double color_b, double  color_g);
//Отрисовка линии
void draw_line(double x_start, double y_start, double z_start,
               double x_end, double y_end, double z_end);

//Параметры точек
void point_settings(double point_size, double color_r, double color_b, double color_g);
//Отрисовка точек
void draw_point(double x_pos, double y_pos, double z_pos);

//Отрисовка текста - подготовка
void draw_text_preparation(double coord_x, double coord_y, double coord_z, double size_coeff);
//Пишем текст
void draw_string_stroke(std::string input);
void draw_string_3d(std::string input, double coord_x, double coord_y, double coord_z, double size_coeff);
//Полный текст
void draw_string(std::string input, double coord_x, double coord_y, double size_coeff);

// Отрисовка осей
void draw_axis();
//Получение порядка числа
int GetMagnitudeOfInteger(double x);
//Отрисовка засечек
void draw_axis_marks();
//Легенда
void draw_legend();
//Сетка
void draw_net();

//Задание размеров изображение после вычисления
void set_bounds();
//Задание размеров изображение после вычисления функции
void set_f_bounds();
//Задание размеров изображения графика для точеченого графика
void initialization_bounds();

//Инициализация для функций
void initialization_function(double i_diapason_x_start, double i_diapason_x_end,
                             double(**i_functions)(double), int i_functions_number);
//Инициализация для поверхностей
void initialization_surface(double i_diapason_x_start, double i_diapason_x_end,
                            double i_diapason_y_start, double i_diapason_y_end,
                            double(**i_surface_functions)(double, double),
                            int i_surface_number);

//Цвет поверхности
//Преобразуем точку y в точку между 0 и 1
//На основе этого вычисляем цвет
void surface_coloring(double axis_start, double axis_end, double current_axis);


//Отрисовка треугольника
void draw_triangle(Point first_point, Point second_point, Point third_point, double t_size);

void draw_surface(int t);


//Начальные настройки графики
void graphic_initialization();
void draw_graphic();

//Отрисовка графика
void draw();

// Отрисовка кадра
void display();

//Изменение размера
void OnRedraw(GLint w, GLint h);
//Замыкание круга
double degree_circling(double x);

//Обработчик клавиатуры. 
void keyboard(unsigned char key, GLint x, GLint y);
// Отрисовка
void OpenGL();


//Инициализация наборов точек
void initialization_p_many(double ***i_points_xyz, int i_point_number_x, int i_axis_number, int i_sets_number);


// ПОЛЬЗОВАТЕЛЬСКИЕ НАСТРОЙКИ

//Оси
void Set_Axis_X_Text(std::string new_text); // Текст по оси x
void Set_Axis_Y_Text(std::string new_text); // Текст по оси y
void Set_Axis_Z_Text(std::string new_text); // Текст по оси z
void Set_Axis_Text_Size(double new_size);   //Размер шрифта текста
void Set_Text_Line_Size(double new_size);   //Размер линии текста
void Set_Marks_Text_Size(double new_size);  //Размер шрифта текста засечек
void Show_Z_Marks(bool boolean);            // Показать засечки по оси Z
void Set_Axis_Point_Size(double new_size);  // Размер точек осей
void Set_Axis_Line_Size(double new_size);   // Размер линий осей
void Set_Marks_Line_Size(double new_size);  // Размер засечек осей 

// График
void surface_mode(bool new_mode);
void Set_Step_X(double new_step);             // Шаг по оси x
void Set_Step_Y(double new_step);             // Шаг по оси y
void Set_Graphic_Line_Size(double new_size);  // Размер линий графика
void Set_Graphic_Point_Size(double new_size); // Размер точек графика
// Поверхность
void Set_Surface_Line_Size(double new_size);
void Set_Surface_Point_Size(double new_size);
void Set_Surface_Width(double new_size);


//Сетка
void Set_Net(bool flag);                    // Нужна сетка?
void Set_Net_Line_Size(double new_size);    // Размер линий сетки
void tenth_mode(bool flag);                 // Режим миллиметра

// Легенда
void Set_Legend_Text(std::string* new_text); // Текст легенды
void Set_Legend_Start_Pos_X(double posit);   // Позиция легенды
void Show_Legend(bool flag);                 // Показать легенду?