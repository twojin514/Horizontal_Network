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
#include <iomanip>
#include <codecvt>
#include <cmath>

struct Point {
    int index;
    std::string name;
    double X, Y, X_SD, Y_SD;
    bool isControl;
};

struct Line {
    int index;
    std::string From, To;
    double Distance, Distance_SD;
};

struct Angle {
    int index;
    std::string Backsight, Occupied, Foresight;
    double Angle_D, Angle_M, Angle_S;
    double Angle_SD;
};

void PrintInformation(std::ofstream& outfile, const std::string& title, const std::string& business_information, const std::string& name, const std::vector<Point>& points_data, const std::vector<Line>& lines_data, const std::vector<Angle>& angles_data);
void ReadPointFile(const std::string& Control_FileName, const std::string& Coordinate_FileName, std::vector<Point>& point_data);
void ReadControlFile(const std::string& Control_FileName, std::vector<Point>& control_data);
void ReadMeasureFile(const std::string& Coordinate_FileName, std::vector<Point>& measure_data);
void ReadLineFile(const std::string& Distance_FileName, std::vector<Line>& lines_data);
void ReadAngleFile(const std::string& Angle_FileName, std::vector<Angle>& angles_data);
bool isNumber(const std::string& s);

std::string removeSpaces(const std::string& str);

void SortFile(std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles);
void PrintSortFile(std::ofstream& outfile, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles);
void ComposeWMatrix(cv::Mat& W, const std::vector<Point>& Control_Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles);
void ComposeAMatrix(cv::Mat& A, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const std::vector<Point>& Control_Points);
void ComposeLMatrix(cv::Mat& L, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const std::vector<Point>& Control_Points, std::vector<double>& predicteDistances, std::vector<double>& predicteAngles);
void CalculateXMatrix(cv::Mat& X, const cv::Mat& A, const cv::Mat& W, const cv::Mat& L);
void CalculateVMatrix(cv::Mat& V, const cv::Mat& A, const cv::Mat& X, const cv::Mat& L);
void CalculateSoMatrix(cv::Mat& So, const cv::Mat& V, const cv::Mat& W, int m, int n);
void PrintIteration(std::ofstream& outfile, int iteration, const std::vector<Point>& Points, const cv::Mat& X, const cv::Mat& V, const cv::Mat& L, const cv::Mat& So, int m, int n);
void UpdateValues(std::vector<Point>& Points, const cv::Mat& X);
void PrintResult(std::ofstream& outfile, int& m, const std::vector<Point>& Points, const std::vector<Point>& Control_Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const cv::Mat& X, const std::vector<Point>& Points_init, const std::vector<double>& predicteDistances, const std::vector<double>& predicteAngles, const cv::Mat& SigmaXX, const cv::Mat& SigmaLL, const cv::Mat& So, const double X2);