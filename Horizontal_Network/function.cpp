#include "function.h"
#define M_PI       3.14159265358979323846   // pi

void ReadPointFile(std::string Control_FileName, std::string Coordinate_FileName, std::vector<Point>& points_data)
{
    std::ifstream file1(Control_FileName);
    std::ifstream file2(Coordinate_FileName);
    std::string point_header;
    std::string point_data;
    std::vector<Point> vetor_points;
    int i = 0;

    if (!file1.is_open())
    {
        printf("Error : Point File not found\n");
        return; // 파일 열기 실패 시 함수 종료
    }

    if (!file2.is_open())
    {
        printf("Error : Point File not found\n");
        return; // 파일 열기 실패 시 함수 종료
    }

    std::getline(file2, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file2, point_data))
    {
        std::istringstream iss(point_data);
        Point point;
        std::string token;
        point.index = ++i;
        std::getline(iss, token, ',');
        point.name = token;
        std::getline(iss, token, ',');
        point.X = std::stod(token);
        std::getline(iss, token, ',');
        point.Y = std::stod(token);
        point.X_SD = NULL;
        point.Y_SD = NULL;
        point.isControl = false;

        vetor_points.push_back(point);
    }

    std::getline(file1, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file1, point_data))
    {
        std::istringstream iss(point_data);
        Point point;
        std::string token;
        point.index = ++i;
        std::getline(iss, token, ',');
        point.name = token;
        std::getline(iss, token, ',');
        point.X = std::stod(token);
        std::getline(iss, token, ',');
        point.Y = token.empty() ? 0 : std::stod(token);
        std::getline(iss, token, ',');
        point.X_SD = token.empty() ? 0 : std::stod(token);
        std::getline(iss, token, ',');
        point.Y_SD = token.empty() ? 0 : std::stod(token);
        point.isControl = true;

        vetor_points.push_back(point);
    }

    


    points_data = vetor_points;

}
void ReadControlFile(std::string Control_FileName, std::vector<Point>& control_data)
{
	std::ifstream file(Control_FileName);
	std::string point_header;
	std::string point_data;
	std::vector<Point> vetor_points;
	int i = 0;

	if (!file.is_open())
	{
		printf("Error : Point File not found\n");
		return; // 파일 열기 실패 시 함수 종료
	}

	std::getline(file, point_header); // 첫 줄은 헤더이므로 읽지 않음

	while (std::getline(file, point_data))
	{
		std::istringstream iss(point_data);
		Point point;
		std::string token;
		point.index = ++i;
		std::getline(iss, token, ',');
		point.name = token;
		std::getline(iss, token, ',');
		point.X = std::stod(token);
		std::getline(iss, token, ',');
		point.Y = token.empty() ? 0 : std::stod(token);
		std::getline(iss, token, ',');
		point.X_SD = token.empty() ? 0 : std::stod(token);
		std::getline(iss, token, ',');
		point.Y_SD = token.empty() ? 0 : std::stod(token);
		point.isControl = true;

		vetor_points.push_back(point);
	}

	control_data = vetor_points;
}
void ReadMeasureFile(std::string Coordinate_FileName, std::vector<Point>& measure_data)
{
	std::ifstream file(Coordinate_FileName);
	std::string point_header;
	std::string point_data;
	std::vector<Point> vetor_points;
	int i = 0;

	if (!file.is_open())
	{
		printf("Error : Point File not found\n");
		return; // 파일 열기 실패 시 함수 종료
	}

	std::getline(file, point_header); // 첫 줄은 헤더이므로 읽지 않음

	while (std::getline(file, point_data))
	{
		std::istringstream iss(point_data);
		Point point;
		std::string token;
		point.index = ++i;
		std::getline(iss, token, ',');
		point.name = token;
		std::getline(iss, token, ',');
		point.X = std::stod(token);
		std::getline(iss, token, ',');
		point.Y = token.empty() ? 0 : std::stod(token);
		std::getline(iss, token, ',');
		point.X_SD = token.empty() ? 0 : std::stod(token);
		std::getline(iss, token, ',');
		point.Y_SD = token.empty() ? 0 : std::stod(token);
		point.isControl = false;

		vetor_points.push_back(point);
	}

	measure_data = vetor_points;
}
void ReadLineFile(std::string Distance_FileName, std::vector<Line>& lines_data)
{
    std::ifstream file(Distance_FileName);
    std::string line_header;
    std::string line_data;
    std::vector<Line> vetor_lines;
    int i = 0;

    if (!file.is_open())
    {
        printf("Error : Line File not found\n");
        return; // 파일 열기 실패 시 함수 종료
    }

    std::getline(file, line_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, line_data))
    {
        std::istringstream iss(line_data);
        Line line;
        std::string token;
        line.index = ++i;
        std::getline(iss, token, ',');
        line.From = token;
        std::getline(iss, token, ',');
        line.To = token;
        std::getline(iss, token, ',');
        line.Distance = std::stod(token);
        std::getline(iss, token, ',');
        line.Distance_SD = token.empty() ? 0 : std::stod(token);

        vetor_lines.push_back(line);
    }

    lines_data = vetor_lines;
}
void ReadAngleFile(std::string Angle_FileName, std::vector<Angle>& angles_data)
{
    std::ifstream file(Angle_FileName);
    std::string angle_header;
    std::string angle_data;
    std::vector<Angle> vetor_angles;
    int i = 0;

    if (!file.is_open())
    {
        printf("Error : Angle File not found\n");
        return; // 파일 열기 실패 시 함수 종료
    }

    std::getline(file, angle_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, angle_data))
    {
        std::istringstream iss(angle_data);
        Angle angle;
        std::string token;
        angle.index = ++i;
        std::getline(iss, token, ',');
        angle.Backsight = token;
        std::getline(iss, token, ',');
        angle.Occupied = token;
        std::getline(iss, token, ',');
        angle.Foresight = token;
        std::getline(iss, token, ',');
        angle.Angle_D = std::stod(token);
        std::getline(iss, token, ',');
        angle.Angle_M = std::stod(token);
        std::getline(iss, token, ',');
        angle.Angle_S = std::stod(token);
        std::getline(iss, token, ',');
        angle.Angle_SD = std::stod(token);

        vetor_angles.push_back(angle);
    }

    angles_data = vetor_angles;
}
void PrintInformation(std::ofstream& outfile, std::string& title, std::string& business_information, std::string& name, std::vector<Point>& points_data, std::vector<Line>& lines_data, std::vector<Angle>& angles_data)
{
    outfile << std::fixed << std::setprecision(3);


    outfile << "**********************************************************************************************\n";
    outfile << "************************************ " << title << " ******************************\n";
    outfile << "**************************** " << business_information << " *************************\n";
    outfile << "**********************************************************************************************\n";
    outfile << "***************************************************** Name : " << name << " **********************\n";

    outfile << "**********************************************************************************************\n\n";
    outfile << "\n*************************************** 1. Input Ddata ***************************************\n";
    outfile << "\n************************************* Point Data **************************************\n";
    outfile << std::setw(15) << "Point_ID" << "\t|"
        << std::setw(25) << "Control(T/F)" << "\t|"
        << std::setw(25) << "Point_Name" << "\t|"
        << std::setw(25) << "X" << "\t|"
        << std::setw(25) << "Y" << "\t|"
        << std::setw(25) << "S_X" << "\t|"
        << std::setw(25) << "S_Y" << "\n\n";


    for (int i = 0; i < points_data.size(); ++i) {
        outfile << std::setw(12) << points_data[i].index << "\t"
            << std::setw(25) << (points_data[i].isControl ? "T" : "F") << "\t"
            << std::setw(25) << points_data[i].name << "\t"
            << std::setw(25) << points_data[i].X << "\t"
            << std::setw(25) << points_data[i].Y << "\t"
            << std::setw(25) << points_data[i].X_SD << "\t"
            << std::setw(25) << points_data[i].Y_SD << "\n\n";

        std::cout << std::setw(12) << (points_data[i].isControl ? "T" : "F") << "\t"
            << std::setw(25) << points_data[i].name << "\t"
            << std::setw(25) << points_data[i].X << "\t"
            << std::setw(25) << points_data[i].Y << "\t"
            << std::setw(25) << (points_data[i].isControl ? " " : " ") << "\t"
            << std::setw(25) << (points_data[i].isControl ? " " : " ") << "\n\n";
    }

    outfile << "\n************************************* Line Data **************************************\n";
    std::cout << "\n************************************* Line Data ****************************************\n";

    outfile << std::setw(10) << "Line_ID" << "\t|"
        << std::setw(21) << "From" << "\t|"
        << std::setw(21) << "To" << "\t|"
        << std::setw(21) << "Distance" << "\t|"
        << std::setw(21) << "Distance_SD" << "\n\n";
    std::cout << std::setw(10) << "Line_ID" << "\t|"
        << std::setw(21) << "From" << "\t|"
        << std::setw(21) << "To" << "\t|"
        << std::setw(21) << "Distance" << "\t|"
        << std::setw(21) << "Distance_SD" << "\n\n";

    for (int i = 0; i < lines_data.size(); ++i) {
        outfile << std::setw(10) << lines_data[i].index << "\t"
            << std::setw(21) << lines_data[i].From << "\t"
            << std::setw(21) << lines_data[i].To << "\t"
            << std::setw(21) << lines_data[i].Distance << "\t"
            << std::setw(21) << lines_data[i].Distance_SD << "\n\n";

        std::cout << std::setw(10) << lines_data[i].index << "\t"
            << std::setw(21) << lines_data[i].From << "\t"
            << std::setw(21) << lines_data[i].To << "\t"
            << std::setw(21) << lines_data[i].Distance << "\t"
            << std::setw(21) << lines_data[i].Distance_SD << "\n\n";
    }

    outfile << "\n************************************* Angle Data **************************************\n";
    std::cout << "\n************************************* Angle Data ****************************************\n";

    outfile << std::setw(10) << "Angle_ID" << "\t|"
        << std::setw(21) << "Backsight" << "\t|"
        << std::setw(21) << "Occupied" << "\t|"
        << std::setw(21) << "Foresight" << "\t|"
        << std::setw(21) << "Angle_D" << "\t|"
        << std::setw(21) << "Angle_S" << "\t|"
        << std::setw(21) << "Angle_M" << "\t|"
        << std::setw(21) << "Angle_SD" << "\n\n";

    for (int i = 0; i < angles_data.size(); ++i) {
        outfile << std::setw(10) << angles_data[i].index << "\t"
            << std::setw(21) << angles_data[i].Backsight << "\t"
            << std::setw(21) << angles_data[i].Occupied << "\t"
            << std::setw(21) << angles_data[i].Foresight << "\t"
            << std::setw(21) << angles_data[i].Angle_D << "\t"
            << std::setw(21) << angles_data[i].Angle_S << "\t"
            << std::setw(21) << angles_data[i].Angle_M << "\t"
            << std::setw(21) << angles_data[i].Angle_SD << "\n\n";
    }
}
void SortFile(std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles)
{
    std::vector<Point> points = Points;
    std::vector<Line> distances = Distances;
    std::vector<Angle> angles = Angles;

    std::vector<std::string> point_names;


    // Collect All Unique Point Names
    for (auto point : points) {
        if (std::find(point_names.begin(), point_names.end(), point.name) == point_names.end()) {
            point_names.push_back(point.name);
        }
    }
    for (auto distance : distances) {
        if (std::find(point_names.begin(), point_names.end(), distance.From) == point_names.end()) {
            point_names.push_back(distance.From);
        }
        if (std::find(point_names.begin(), point_names.end(), distance.To) == point_names.end()) {
            point_names.push_back(distance.To);
        }
    }
    for (auto angle : angles) {
        if (std::find(point_names.begin(), point_names.end(), angle.Backsight) == point_names.end()) {
            point_names.push_back(angle.Backsight);
        }
        if (std::find(point_names.begin(), point_names.end(), angle.Occupied) == point_names.end()) {
            point_names.push_back(angle.Occupied);
        }
        if (std::find(point_names.begin(), point_names.end(), angle.Foresight) == point_names.end()) {
            point_names.push_back(angle.Foresight);
        }
    }

    // Sort The Point Names in Alphabetical Order
    for (auto it = point_names.begin(); it != point_names.end(); ) {
        if (*it == "") {
            it = point_names.erase(it);
        }
        else {
            ++it;
        }
    }
    std::sort(point_names.begin(), point_names.end());

    // Save Sorting Result
    Points = points;
    Distances = distances;
    Angles = angles;
}

void PrintSortFile(std::ofstream& outfile, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles)
{
    std::vector<Point> points = Points;
    std::vector<Line> distances = Distances;
    std::vector<Angle> angles = Angles;
    outfile << std::fixed << std::setprecision(3);

    outfile << "\n*************************************** 2. Sort station names and order logic ***************************************\n";
    outfile << "\n************************************* Point Data **************************************\n";
    outfile << std::setw(15) << "Point_Order" << "\t|"
        << std::setw(25) << "Point_Name" << "\t|"
        << std::setw(25) << "X" << "\t|"
        << std::setw(25) << "Y" << "\t|"
        << std::setw(25) << "S_X" << "\t|"
        << std::setw(25) << "S_Y" << "\n\n";

    for (int i = 0; i < points.size(); ++i) {
        outfile << std::setw(12) << points[i].index << "\t"
            << std::setw(25) << points[i].name << "\t"
            << std::setw(25) << points[i].X << "\t"
            << std::setw(25) << points[i].Y << "\t"
            << std::setw(25) << points[i].X_SD << "\t"
            << std::setw(25) << points[i].Y_SD << "\n\n";
    }
}

void ComposeWMatrix(cv::Mat& W, int& m, int& n, std::vector<Point>& Control_Points, std::vector<Line>& Distances, std::vector<Angle>& Angles)
{
    std::vector<Point> points = Control_Points;
    std::vector<Line> distances = Distances;
    std::vector<Angle> angles = Angles;

    // Distances
    for (int i = 0; i < Distances.size(); ++i) {
        W.at<double>(i, i) = 1 / (Distances[i].Distance_SD * Distances[i].Distance_SD);
    }

    // Angles
    for (int i = 0, j = Distances.size(); i < Angles.size(); i++, j++) 
		{
			W.at<double>(j, j) = 1 / (Angles[i].Angle_SD/3600 * Angles[i].Angle_SD/3600);
		}

	// Control Points
	for (int i = 0, j = Distances.size() + Angles.size(); i < Control_Points.size(); i++, j++)
	{
		W.at<double>(j, j) = 1 / (Control_Points[i].X_SD * Control_Points[i].X_SD);
		W.at<double>(j + 1, j + 1) = 1 / (Control_Points[i].Y_SD * Control_Points[i].Y_SD);
		j++;
    }
}

void ComposeAMatrix(cv::Mat& A, int& m, int& n, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, std::vector<Point>& Control_Points) {
    std::vector<Point> points = Points;
    std::vector<Line> distances = Distances;
    std::vector<Angle> angles = Angles;
    int q, p, e;

    // Distances
    for (int i = 0; i < Distances.size(); ++i) {
        int fromIdx = -1;
        int toIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Distances[i].From == Points[j].name) {
                fromIdx = j;
            }
            if (Distances[i].To == Points[j].name) {
                toIdx = j;
            }
        }

        if (fromIdx == -1 || toIdx == -1) {
            std::cerr << "Error: Point not found for distance " << Distances[i].From << " to " << Distances[i].To << std::endl;
            continue;
        }

        double deltaX = Points[toIdx].X - Points[fromIdx].X;
        double deltaY = Points[toIdx].Y - Points[fromIdx].Y;
        double distance = sqrt(deltaX * deltaX + deltaY * deltaY);

        A.at<double>(i, 2 * fromIdx) = deltaX / distance;
        A.at<double>(i, 2 * fromIdx + 1) = deltaY / distance;
        A.at<double>(i, 2 * toIdx) = -deltaX / distance;
        A.at<double>(i, 2 * toIdx + 1) = -deltaY / distance;
    }

    // Angles
    for (int i = 0; i < Angles.size(); ++i) {
        int backIdx = -1;
        int occIdx = -1;
        int foreIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Angles[i].Backsight == Points[j].name) {
                backIdx = j;
            }
            if (Angles[i].Occupied == Points[j].name) {
                occIdx = j;
            }
            if (Angles[i].Foresight == Points[j].name) {
                foreIdx = j;
            }
        }

        if (backIdx == -1 || occIdx == -1 || foreIdx == -1) {
            std::cerr << "Error: Point not found for angle " << Angles[i].Backsight << ", " << Angles[i].Occupied << ", " << Angles[i].Foresight << std::endl;
            continue;
        }

        double deltaXBackOcc = Points[backIdx].X - Points[occIdx].X;
        double deltaYBackOcc = Points[backIdx].Y - Points[occIdx].Y;
        double deltaXForeOcc = Points[foreIdx].X - Points[occIdx].X;
        double deltaYForeOcc = Points[foreIdx].Y - Points[occIdx].Y;

        double denomBack = pow(deltaXBackOcc, 2) + pow(deltaYBackOcc, 2);
        double denomFore = pow(deltaXForeOcc, 2) + pow(deltaYForeOcc, 2);

        A.at<double>(i + Distances.size(), 2 * backIdx) = -deltaYBackOcc / denomBack;
        A.at<double>(i + Distances.size(), 2 * backIdx + 1) = deltaXBackOcc / denomBack;
        A.at<double>(i + Distances.size(), 2 * occIdx) = (deltaXBackOcc / denomBack - deltaXForeOcc / denomFore);
        A.at<double>(i + Distances.size(), 2 * occIdx + 1) = (deltaYBackOcc / denomBack - deltaYForeOcc / denomFore);
        A.at<double>(i + Distances.size(), 2 * foreIdx) = deltaYForeOcc / denomFore;
        A.at<double>(i + Distances.size(), 2 * foreIdx + 1) = -deltaXForeOcc / denomFore;
    }

    // Control Points
    for (int i = 0; i < Control_Points.size(); ++i) {
        int ctrlIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Control_Points[i].name == Points[j].name) {
                ctrlIdx = j;
                break;
            }
        }

        if (ctrlIdx == -1) {
            std::cerr << "Error: Control point not found " << Control_Points[i].name << std::endl;
            continue;
        }

        int rowIdx = Distances.size() + Angles.size() + 2 * i;
        A.at<double>(rowIdx, 2 * ctrlIdx) = 1;
        A.at<double>(rowIdx + 1, 2 * ctrlIdx + 1) = 1;
    }
}


void ComposeLMatrix(cv::Mat& L, int& m, std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, std::vector<Point>& Control_Points, std::vector<double>& predicteDistances,std::vector<double>& predicteAngles) {
    std::vector<Point> points = Points;
    std::vector<Line> distances = Distances;
    std::vector<Angle> angles = Angles;

    // Distances
    for (int i = 0; i < Distances.size(); ++i) {
        double observedDistance = Distances[i].Distance;
        double predictedDistance = 0.0;

        int fromIdx = -1;
        int toIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Distances[i].From == Points[j].name) {
                fromIdx = j;
            }
            if (Distances[i].To == Points[j].name) {
                toIdx = j;
            }
        }

        if (fromIdx != -1 && toIdx != -1) {
            double deltaX = Points[toIdx].X - Points[fromIdx].X;
            double deltaY = Points[toIdx].Y - Points[fromIdx].Y;
            predictedDistance = sqrt(deltaX * deltaX + deltaY * deltaY);
        }
        L.at<double>(i, 0) = observedDistance - predictedDistance;
		predicteDistances.push_back(predictedDistance);

    }

    // Angles
    for (int i = 0; i < Angles.size(); ++i) {
        double observedAngle = Angles[i].Angle_D + Angles[i].Angle_M / 60.0 + Angles[i].Angle_S / 3600.0;
        double predictedAngle = 0.0;

        int backIdx = -1;
        int occIdx = -1;
        int foreIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Angles[i].Backsight == Points[j].name) {
                backIdx = j;
            }
            if (Angles[i].Occupied == Points[j].name) {
                occIdx = j;
            }
            if (Angles[i].Foresight == Points[j].name) {
                foreIdx = j;
            }
        }

        double angle;

		if (backIdx != -1 && occIdx != -1 && foreIdx != -1) {
			double deltabf = sqrt(pow(Points[backIdx].X - Points[foreIdx].X, 2) + pow(Points[backIdx].Y - Points[foreIdx].Y, 2));
			double deltabo = sqrt(pow(Points[backIdx].X - Points[occIdx].X, 2) + pow(Points[backIdx].Y - Points[occIdx].Y, 2));
			double deltaof = sqrt(pow(Points[occIdx].X - Points[foreIdx].X, 2) + pow(Points[occIdx].Y - Points[foreIdx].Y, 2));

            double temp = (pow(deltabo, 2) + pow(deltaof, 2) - pow(deltabf, 2)) / (2 * deltabo * deltaof);

            angle = acos(temp) * (180 / M_PI);


            if (observedAngle > 180) {
				angle = 360 - angle;
            }

            angle = angle / (180 / M_PI);

		}
        L.at<double>(i + Distances.size(), 0) = observedAngle/(180/M_PI) - angle;

		double angle_DD = angle * (180 / M_PI);
		predicteAngles.push_back(angle_DD);

    }

    // Control Points
    for (int i = 0; i < Control_Points.size(); ++i) {
        int ctrlIdx = -1;

        for (int j = 0; j < Points.size(); ++j) {
            if (Control_Points[i].name == Points[j].name) {
                ctrlIdx = j;
                break;
            }
        }

        if (ctrlIdx != -1) {
            int rowIdx = Distances.size() + Angles.size() + 2 * i;
            L.at<double>(rowIdx, 0) = Points[ctrlIdx].X - Control_Points[i].X;
            L.at<double>(rowIdx + 1, 0) = Points[ctrlIdx].Y - Control_Points[i].Y;
        }
    }
}

void CalculateXMatrix(cv::Mat& X, cv::Mat& A, cv::Mat& W, cv::Mat& L, int& m, int& n)
{
	cv::Mat AtWA = A.t() * W * A;
	cv::Mat AtWL = A.t() * W * L;
	cv::Mat AtWA_inv = AtWA.inv();
	X = AtWA.inv() * AtWL;
}

void CalculateVMatrix(cv::Mat& V, cv::Mat& A, cv::Mat& X, cv::Mat& L, int& m, int& n)
{
	V = A * X - L;
}

void CalculateSoMatrix(cv::Mat& So, cv::Mat& V, cv::Mat& W, cv::Mat& X, int& m, int& n)
{
	So = (V.t() * W * V)/(m-n);

}

void PrintIteration(std::ofstream& outfile, int& iteration, std::vector<Point>& Points, cv::Mat& X, cv::Mat& V, cv::Mat& L, cv::Mat& So, int& m, int& n)
{
	outfile << std::fixed << std::setprecision(3);

	outfile << "\n*************************************** Iteration " << iteration +1 << " ***************************************\n\n";
	outfile << std::setw(12) << "index" << "\t|"
		<< std::setw(25) << "Point" << "\t|"
		<< std::setw(25) << "before_X" << "\t|"
		<< std::setw(25) << "before_Y" << "\t|"
		<< std::setw(25) << "dX" << "\t|"
		<< std::setw(25) << "dY" << "\t|"
		<< std::setw(25) << "after_X" << "\t|"
		<< std::setw(25) << "after_Y" << "\n\n";
    
    for (int i = 0; i < n / 2; ++i) {
		outfile << std::setw(12) << i << "\t"
			<< std::setw(25) << Points[i].name << "\t"
			<< std::setw(25) << Points[i].X << "\t"
			<< std::setw(25) << Points[i].Y << "\t"
			<< std::setw(25) << X.at<double>(2 * i, 0) << "\t"
			<< std::setw(25) << X.at<double>(2 * i + 1, 0) << "\t"
			<< std::setw(25) << Points[i].X + X.at<double>(2 * i, 0) << "\t"
			<< std::setw(25) << Points[i].Y + X.at<double>(2 * i + 1, 0) << "\n\n";
    }



	outfile << "\n************************************* Variance-Covariance Matrix **************************************\n";
	outfile << std::setw(10) << "So : " << sqrt(So.at<double>(0, 0)) << "\n\n";
	
}

void UpdateValues(std::vector<Point>& Points, cv::Mat& X, int& n)
{
	for (int i = 0; i < n/2; ++i) {
		Points[i].X += X.at<double>(2*i, 0);
		Points[i].Y += X.at<double>(2*i + 1, 0);
	}


}

void PrintResult(std::ofstream& outfile, std::vector<Point>& Points, std::vector<Point>& Control_Points, std::vector<Line>& Distances, std::vector<Angle>& Angles, cv::Mat& X, std::vector<Point>& Points_init, int& n, std::vector<double>& predicteDistances, std::vector<double>& predicteAngles)
{
	outfile << std::fixed << std::setprecision(3);

	outfile << "\n*************************************** End ***************************************\n\n";

	outfile << "\n************************************* Point Data **************************************\n";
	outfile << std::setw(15) << "Point_ID" << "\t|"
		<< std::setw(25) << "Point_Name" << "\t|"
		<< std::setw(25) << "X_init" << "\t|"
		<< std::setw(25) << "Y_init" << "\t|"
		<< std::setw(25) << "X_final" << "\t|"
		<< std::setw(25) << "Y_final" << "\t|"
		<< std::setw(25) << "SD_X" << "\t|"
		<< std::setw(25) << "SD_Y" << "\n\n";

    for (int i = 0; i < Points.size(); ++i) {
		outfile << std::setw(12) << Points[i].index << "\t"
			<< std::setw(25) << Points[i].name << "\t"
			<< std::setw(25) << Points_init[i].X << "\t"
			<< std::setw(25) << Points_init[i].Y << "\t"
			<< std::setw(25) << Points[i].X << "\t"
			<< std::setw(25) << Points[i].Y << "\t"
			<< std::setw(25) << Points[i].X_SD << "\t"
			<< std::setw(25) << Points[i].Y_SD << "\n\n";
    }

	outfile << "\n************************************* Line Data **************************************\n";
	outfile << std::setw(15) << "Line_ID" << "\t|"
		<< std::setw(25) << "From" << "\t|"
		<< std::setw(25) << "To" << "\t|"
		<< std::setw(25) << "Distance_obs" << "\t|"
		<< std::setw(25) << "Distance_cal" << "\t|"
		<< std::setw(25) << "V" << "\n\n";

    for (int i = 0; i < Distances.size(); ++i) {
        outfile << std::setw(15) << Distances[i].index << "\t"
            << std::setw(25) << Distances[i].From << "\t"
            << std::setw(25) << Distances[i].To << "\t"
            << std::setw(25) << Distances[i].Distance << "\t"
            << std::setw(25) << predicteDistances[i] << "\t"
            << std::setw(25) << Distances[i].Distance - predicteDistances[i] << "\n\n";
    }

	outfile << "\n************************************* Angle Data **************************************\n";
	outfile << std::setw(15) << "Angle_ID" << "\t|"
		<< std::setw(25) << "Backsight" << "\t|"
		<< std::setw(25) << "Occupied" << "\t|"
		<< std::setw(25) << "Foresight" << "\t|"
		<< std::setw(25) << "Angle_obs" << "\t|"
		<< std::setw(25) << "Angle_cal" << "\t|"
		<< std::setw(25) << "V" << "\n\n";

	for (int i = 0; i < Angles.size(); ++i) {
		outfile << std::setw(15) << Angles[i].index << "\t"
			<< std::setw(25) << Angles[i].Backsight << "\t"
			<< std::setw(25) << Angles[i].Occupied << "\t"
			<< std::setw(25) << Angles[i].Foresight << "\t"
			<< std::setw(25) << Angles[i].Angle_D + Angles[i].Angle_M / 60.0 + Angles[i].Angle_S / 3600.0 << "\t"
			<< std::setw(25) << predicteAngles[i] << "\t"
			<< std::setw(25) << (Angles[i].Angle_D + Angles[i].Angle_M / 60.0 + Angles[i].Angle_S / 3600.0) - predicteAngles[i] << "\n\n";
	}




}






