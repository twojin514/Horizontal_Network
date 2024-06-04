#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <locale>
#include <cctype>
#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <cctype>
#include <locale>
#include <iomanip>
#include <codecvt>

struct Point
{
	int index;
	std::string name;
	double X, Y, X_SD, Y_SD;
	bool isControl;
};

struct Line
{
	int index;
	std::string From, To;
	double Distance, Distance_SD;
};

struct Angle
{
	int index;
	std::string Backsight, Occupied, Foresight;
	double Angle_D, Angle_M, Angle_S;
	double Angle_Deg, Angle_SD;
};

void PrintInformation(std::ofstream& outfile, std::string& title, std::string& business_information, std::string& name, std::vector<Point>& points_data, std::vector<Line>& lines_data, std::vector<Angle>& angles_data);
void ReadPointFile(std::string Control_FileName, std::string Coordinate_FileName, std::vector<Point>& point_data);
void ReadControlFile(std::string Control_FileName, std::vector<Point>& control_data);
void ReadMeasureFile(std::string Coordinate_FileName, std::vector<Point>& measure_data);
void ReadLineFile(std::string Distance_FileName, std::vector<Line>& lines_data);

void ReadAngleFile(std::string Angle_FileName, std::vector<Angle>& angles_data);

void PrintInformation(std::ofstream& outfile, std::string& title, std::string& business_information, std::string& name, std::vector<Point>& points_data, std::vector<Line>& lines_data, std::vector<Angle>& angles_data);



void SortFile(std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles);
void PrintSortFile(std::ofstream& outfile, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles);

void ComposeWMatrix(cv::Mat& W, int& m, int& n, std::vector<Point>& Control_Points, std::vector<Line>& Distances, std::vector<Angle>& Angles);
void ComposeAMatrix(cv::Mat& A, int& m, int& n, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, std::vector<Point>& Control_Points);
void ComposeLMatrix(cv::Mat& L, int& m, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, std::vector<Point>& Control_Points, std::vector<double>& predicteDistances, std::vector<double>& predicteAngles);
	void CalculateXMatrix(cv::Mat& X, cv::Mat& A, cv::Mat& W, cv::Mat& L, int& m, int& n);
void CalculateVMatrix(cv::Mat& V, cv::Mat& A, cv::Mat& X, cv::Mat& L, int& m, int& n);
void CalculateSoMatrix(cv::Mat& So, cv::Mat& V, cv::Mat& W, cv::Mat& X, int& m, int& n);
void PrintIteration(std::ofstream& outfile, int& iteration, std::vector<Point>& Points, cv::Mat& X, cv::Mat& V, cv::Mat& L, cv::Mat& So, int& m, int& n);
void UpdateValues(std::vector<Point>& Points, cv::Mat& X, int& n);
void PrintResult(std::ofstream& outfile, std::vector<Point>& Points, std::vector<Point>& Control_Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, cv::Mat& X, std::vector<Point>& Points_init, int& n, std::vector<double>& predicteDistances, std::vector<double>& predicteAngles);
