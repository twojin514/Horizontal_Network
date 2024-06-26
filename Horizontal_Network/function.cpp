#include "function.h"
#define M_PI 3.14159265358979323846

void ReadPointFile(const std::string& Control_FileName, const std::string& Coordinate_FileName, std::vector<Point>& points_data) {
    std::ifstream file1(Control_FileName);
    std::ifstream file2(Coordinate_FileName);
    std::string point_header;
    std::string point_data;
    std::vector<Point> vetor_points;
    int i = 0;

    if (!file1.is_open()) {
        std::cerr << "Error : Control File not found\n";
        return;
    }

    if (!file2.is_open()) {
        std::cerr << "Error : Coordinate File not found\n";
        return;
    }

    std::getline(file2, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file2, point_data)) {
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
        point.X_SD = 0;
        point.Y_SD = 0;
        point.isControl = false;

        vetor_points.push_back(point);
    }

    std::getline(file1, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file1, point_data)) {
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

void ReadControlFile(const std::string& Control_FileName, std::vector<Point>& control_data) {
    std::ifstream file(Control_FileName);
    std::string point_header;
    std::string point_data;
    std::vector<Point> vetor_points;
    int i = 0;

    if (!file.is_open()) {
        std::cerr << "Error : Control File not found\n";
        return;
    }

    std::getline(file, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, point_data)) {
        std::istringstream iss(point_data);
        Point point;
        std::string token;
        point.index = ++i;
        std::getline(iss, token, ',');
        point.name = removeSpaces(token);
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

void ReadMeasureFile(const std::string& Coordinate_FileName, std::vector<Point>& measure_data) {
    std::ifstream file(Coordinate_FileName);
    std::string point_header;
    std::string point_data;
    std::vector<Point> vetor_points;
    int i = 0;

    if (!file.is_open()) {
        std::cerr << "Error : Coordinate File not found\n";
        return;
    }

    std::getline(file, point_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, point_data)) {
        std::istringstream iss(point_data);
        Point point;
        std::string token;
        point.index = ++i;
        std::getline(iss, token, ',');
        point.name = removeSpaces(token);
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

//void ReadLineFile(const std::string& Distance_FileName, std::vector<Line>& lines_data) {
//    std::ifstream file(Distance_FileName);
//    std::string line_header;
//    std::string line_data;
//    std::vector<Line> vetor_lines;
//    int i = 0;
//
//    if (!file.is_open()) {
//        std::cerr << "Error : Line File not found\n";
//        return;
//    }
//
//    std::getline(file, line_header); // 첫 줄은 헤더이므로 읽지 않음
//
//    while (std::getline(file, line_data)) {
//        std::istringstream iss(line_data);
//        Line line;
//        std::string token;
//        line.index = ++i;
//        std::getline(iss, token, ',');
//        line.From = token;
//        std::getline(iss, token, ',');
//        line.To = token;
//        std::getline(iss, token, ',');
//        line.Distance = std::stod(token);
//        std::getline(iss, token, ',');
//        line.Distance_SD = token.empty() ? 0 : std::stod(token);
//
//        vetor_lines.push_back(line);
//    }
//
//    lines_data = vetor_lines;
//}
void ReadLineFile(const std::string& Distance_FileName, std::vector<Line>& lines_data) {
    std::ifstream file(Distance_FileName);
    std::string line_header;
    std::string line_data;
    std::vector<Line> vetor_lines;
    int i = 0;

    if (!file.is_open()) {
        std::cerr << "Error : Line File not found\n";
        return;
    }

    std::getline(file, line_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, line_data)) {
        std::istringstream iss(line_data);
        Line line;
        std::string token;
        line.index = ++i;

        // 공백 제거 후 From 필드에 저장
        std::getline(iss, token, ',');
        line.From = removeSpaces(token);

        // 공백 제거 후 To 필드에 저장
        std::getline(iss, token, ',');
        line.To = removeSpaces(token);

        // Distance 필드에 저장
        std::getline(iss, token, ',');
        line.Distance = std::stod(token);

        // Distance_SD 필드에 저장
        std::getline(iss, token, ',');
        line.Distance_SD = token.empty() ? 0 : std::stod(token);

        vetor_lines.push_back(line);
    }

    lines_data = vetor_lines;
}
void ReadAngleFile(const std::string& Angle_FileName, std::vector<Angle>& angles_data) {
    std::ifstream file(Angle_FileName);
    std::string angle_header;
    std::string angle_data;
    std::vector<Angle> vetor_angles;
    int i = 0;

    if (!file.is_open()) {
        std::cerr << "Error : Angle File not found\n";
        return;
    }

    std::getline(file, angle_header); // 첫 줄은 헤더이므로 읽지 않음

    while (std::getline(file, angle_data)) {
        std::istringstream iss(angle_data);
        Angle angle;
        std::string token;
        angle.index = ++i;
        std::getline(iss, token, ',');
        angle.Backsight = removeSpaces(token);
        std::getline(iss, token, ',');
        angle.Occupied = removeSpaces(token);
        std::getline(iss, token, ',');
        angle.Foresight = removeSpaces(token);
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

void PrintInformation(std::ofstream& outfile, const std::string& title, const std::string& business_information, const std::string& name, const std::vector<Point>& points_data, const std::vector<Line>& lines_data, const std::vector<Angle>& angles_data) {
    outfile << std::fixed << std::setprecision(3);

    outfile << "**********************************************************************************************\n";
    outfile << "*************************************** " << title << " **********************************\n";
    outfile << "********************************* " << business_information << " ***************************\n";
    outfile << "**********************************************************************************************\n";
    outfile << "*********************************************************** Name : " << name << " **********************\n";
    outfile << "**********************************************************************************************\n\n";

    outfile << "\n*************************************** 1. Input Data ***************************************\n";
    outfile << "\n**************************************** Point Data *****************************************\n";
    outfile << std::setw(15) << "Point_ID" << "\t|"
        << std::setw(25) << "Control(T/F)" << "\t|"
        << std::setw(25) << "Point_Name" << "\t|"
        << std::setw(25) << "X" << "\t|"
        << std::setw(25) << "Y" << "\t|"
        << std::setw(25) << "S_X" << "\t|"
        << std::setw(25) << "S_Y" << "\n\n";

    for (const auto& point : points_data) {
        outfile << std::setw(12) << point.index << "\t"
            << std::setw(25) << (point.isControl ? "T" : "F") << "\t"
            << std::setw(25) << point.name << "\t"
            << std::setw(25) << point.X << "\t"
            << std::setw(25) << point.Y << "\t"
            << std::setw(25) << point.X_SD << "\t"
            << std::setw(25) << point.Y_SD << "\n\n";
    }

    outfile << "\n***************************************** Line Data *****************************************\n";

    outfile << std::setw(10) << "Line_ID" << "\t|"
        << std::setw(21) << "From" << "\t|"
        << std::setw(21) << "To" << "\t|"
        << std::setw(21) << "Distance" << "\t|"
        << std::setw(21) << "S_Distance" << "\n\n";

    for (const auto& line : lines_data) {
        outfile << std::setw(10) << line.index << "\t"
            << std::setw(21) << line.From << "\t"
            << std::setw(21) << line.To << "\t"
            << std::setw(21) << line.Distance << "\t"
            << std::setw(21) << line.Distance_SD << "\n\n";
    }

    outfile << "\n**************************************** Angle Data *****************************************\n";

    outfile << std::setw(10) << "Angle_ID" << "\t|"
        << std::setw(21) << "Backsight" << "\t|"
        << std::setw(21) << "Occupied" << "\t|"
        << std::setw(21) << "Foresight" << "\t|"
        << std::setw(21) << "Angle_D" << "\t|"
        << std::setw(21) << "Angle_M" << "\t|"
        << std::setw(21) << "Angle_S" << "\t|"
        << std::setw(21) << "S_Angle" << "\n\n";

    for (const auto& angle : angles_data) {
        outfile << std::setw(10) << angle.index << "\t"
            << std::setw(21) << angle.Backsight << "\t"
            << std::setw(21) << angle.Occupied << "\t"
            << std::setw(21) << angle.Foresight << "\t"
            << std::setw(21) << angle.Angle_D << "\t"
            << std::setw(21) << angle.Angle_M << "\t"
            << std::setw(21) << angle.Angle_S << "\t"
            << std::setw(21) << angle.Angle_SD << "\n\n";
    }
}

// 공백 제거 함수
std::string removeSpaces(const std::string& str) {
    std::string result = str;
    result.erase(std::remove(result.begin(), result.end(), ' '), result.end());
    return result;
}

//void SortFile(std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles) {
//    std::vector<std::string> point_names;
//
//    for (const auto& point : Points) {
//        if (std::find(point_names.begin(), point_names.end(), point.name) == point_names.end()) {
//            point_names.push_back(point.name);
//        }
//    }
//    for (const auto& distance : Distances) {
//        if (std::find(point_names.begin(), point_names.end(), distance.From) == point_names.end()) {
//            point_names.push_back(distance.From);
//        }
//        if (std::find(point_names.begin(), point_names.end(), distance.To) == point_names.end()) {
//            point_names.push_back(distance.To);
//        }
//    }
//    for (const auto& angle : Angles) {
//        if (std::find(point_names.begin(), point_names.end(), angle.Backsight) == point_names.end()) {
//            point_names.push_back(angle.Backsight);
//        }
//        if (std::find(point_names.begin(), point_names.end(), angle.Occupied) == point_names.end()) {
//            point_names.push_back(angle.Occupied);
//        }
//        if (std::find(point_names.begin(), point_names.end(), angle.Foresight) == point_names.end()) {
//            point_names.push_back(angle.Foresight);
//        }
//    }
//
//    // Remove empty names and sort
//    point_names.erase(std::remove(point_names.begin(), point_names.end(), ""), point_names.end());
//    std::sort(point_names.begin(), point_names.end());
//
//    // Save Sorting Result
//    Points = Points;
//    Distances = Distances;
//    Angles = Angles;
//}


void SortFile(std::vector<Point>& Points, std::vector<Line>& Distances, std::vector<Angle>& Angles) {
    std::vector<std::string> point_names;

    for (const auto& point : Points) {
        if (std::find(point_names.begin(), point_names.end(), point.name) == point_names.end()) {
            point_names.push_back(point.name);
        }
    }
    for (const auto& distance : Distances) {
        if (std::find(point_names.begin(), point_names.end(), distance.From) == point_names.end()) {
            point_names.push_back(distance.From);
        }
        if (std::find(point_names.begin(), point_names.end(), distance.To) == point_names.end()) {
            point_names.push_back(distance.To);
        }
    }
    for (const auto& angle : Angles) {
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

    // Remove empty names and sort
    point_names.erase(std::remove(point_names.begin(), point_names.end(), ""), point_names.end());
    std::sort(point_names.begin(), point_names.end());

	if (Points.size() != point_names.size()) {
		std::cerr << "Error: Point names do not match" << std::endl;
		return;
	}

    // Sort Points based on point_names
    //std::sort(Points.begin(), Points.end(), [&point_names](const Point& a, const Point& b) {
    //    return std::find(point_names.begin(), point_names.end(), a.name) < std::find(point_names.begin(), point_names.end(), b.name);
    //    });
}

void PrintSortFile(std::ofstream& outfile, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles) {
    outfile << std::fixed << std::setprecision(3);

    outfile << "\n************************* 2. Sort station names and order logic *****************************\n";
    outfile << "\n**************************************** Point Data *****************************************\n";
    outfile << std::setw(15) << "Point_Order" << "\t|"
        << std::setw(25) << "Point_Name" << "\t|"
        << std::setw(25) << "X" << "\t|"
        << std::setw(25) << "Y" << "\t|"
        << std::setw(25) << "S_X" << "\t|"
        << std::setw(25) << "S_Y" << "\n\n";

    for (const auto& point : Points) {
        outfile << std::setw(12) << point.index << "\t"
            << std::setw(25) << point.name << "\t"
            << std::setw(25) << point.X << "\t"
            << std::setw(25) << point.Y << "\t"
            << std::setw(25) << point.X_SD << "\t"
            << std::setw(25) << point.Y_SD << "\n\n";
    }
}

void ComposeWMatrix(cv::Mat& W, const std::vector<Point>& Control_Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles) {
    int m = Distances.size() + Angles.size() + Control_Points.size() * 2;
    W = cv::Mat::zeros(m, m, CV_64F);

    // Distances
    for (int i = 0; i < Distances.size(); ++i) {
        if (Distances[i].Distance_SD > 0) {
            W.at<double>(i, i) = 1 / (Distances[i].Distance_SD * Distances[i].Distance_SD);
        }
        else {
            std::cerr << "Invalid standard deviation for distance at index " << i << std::endl;
        }
    }

    // Angles
    for (int i = 0, j = Distances.size(); i < Angles.size(); i++, j++) {
        double angleSD_rad = (Angles[i].Angle_SD / 3600.0) * (M_PI / 180.0);
        if (angleSD_rad > 0) {
            W.at<double>(j, j) = 1 / (angleSD_rad * angleSD_rad);
        }
        else {
            std::cerr << "Invalid standard deviation for angle at index " << i << std::endl;
        }
    }

    // Control Points
    for (int i = 0, j = Distances.size() + Angles.size(); i < Control_Points.size(); i++, j++) {
        if (Control_Points[i].X_SD > 0) {
            W.at<double>(j, j) = 1 / (Control_Points[i].X_SD * Control_Points[i].X_SD);
        }
        else {
            std::cerr << "Invalid standard deviation for control point X at index " << i << std::endl;
        }
        if (Control_Points[i].Y_SD > 0) {
            W.at<double>(j + 1, j + 1) = 1 / (Control_Points[i].Y_SD * Control_Points[i].Y_SD);
        }
        else {
            std::cerr << "Invalid standard deviation for control point Y at index " << i << std::endl;
        }
        j++;
    }
}

void ComposeAMatrix(cv::Mat& A, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const std::vector<Point>& Control_Points) {
    int m = Distances.size() + Angles.size() + Control_Points.size() * 2;
    int n = Points.size() * 2;
    A = cv::Mat::zeros(m, n, CV_64F);

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

        if (distance > 0) {
            A.at<double>(i, 2 * fromIdx) = -deltaX / distance;
            A.at<double>(i, 2 * fromIdx + 1) = -deltaY / distance;
            A.at<double>(i, 2 * toIdx) = deltaX / distance;
            A.at<double>(i, 2 * toIdx + 1) = deltaY / distance;
        }
        else {
            std::cerr << "Error: Zero distance for points " << Distances[i].From << " and " << Distances[i].To << std::endl;
        }
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



        double deltaXBackOcc = Points[occIdx].X - Points[backIdx].X;
        double deltaYBackOcc = Points[occIdx].Y - Points[backIdx].Y;
        double deltaXForeOcc = -Points[foreIdx].X + Points[occIdx].X;
        double deltaYForeOcc = -Points[foreIdx].Y + Points[occIdx].Y;

        double dBack = deltaXBackOcc * deltaXBackOcc + deltaYBackOcc * deltaYBackOcc;
        double dFore = deltaXForeOcc * deltaXForeOcc + deltaYForeOcc * deltaYForeOcc;

        if (dBack > 0 && dFore > 0) {
            A.at<double>(i + Distances.size(), 2 * backIdx) = deltaYBackOcc / dBack;
            A.at<double>(i + Distances.size(), 2 * backIdx + 1) = -deltaXBackOcc / dBack;

            A.at<double>(i + Distances.size(), 2 * occIdx) = (-deltaYBackOcc / dBack + deltaYForeOcc / dFore);
            A.at<double>(i + Distances.size(), 2 * occIdx + 1) = (deltaXBackOcc / dBack - deltaXForeOcc / dFore);

            A.at<double>(i + Distances.size(), 2 * foreIdx) = -deltaYForeOcc / dFore;
            A.at<double>(i + Distances.size(), 2 * foreIdx + 1) = deltaXForeOcc / dFore;
        }
        else {
            std::cerr << "Error: Zero denominator for angle calculation at index " << i << std::endl;
        }


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

void ComposeLMatrix(cv::Mat& L, const std::vector<Point>& Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const std::vector<Point>& Control_Points, std::vector<double>& predicteDistances, std::vector<double>& predicteAngles) {
    int m = Distances.size() + Angles.size() + Control_Points.size() * 2;
    L = cv::Mat::zeros(m, 1, CV_64F);

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

            angle = angle * (M_PI / 180.0);
        }
        L.at<double>(i + Distances.size(), 0) = observedAngle * (M_PI / 180.0) - angle;

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
            L.at<double>(rowIdx, 0) = Control_Points[i].X - Points[ctrlIdx].X;
            L.at<double>(rowIdx + 1, 0) = Control_Points[i].Y - Points[ctrlIdx].Y;
        }
    }
}

void CalculateXMatrix(cv::Mat& X, const cv::Mat& A, const cv::Mat& W, const cv::Mat& L) {
    cv::Mat AtWA = A.t() * W * A;
    cv::Mat AtWL = A.t() * W * L;
    cv::Mat AtWA_inv = AtWA.inv();
    X = AtWA_inv * AtWL;
}

void CalculateVMatrix(cv::Mat& V, const cv::Mat& A, const cv::Mat& X, const cv::Mat& L) {
    V = A * X - L;
}

void CalculateSoMatrix(cv::Mat& So, const cv::Mat& V, const cv::Mat& W, int m, int n) {
    So = (V.t() * W * V) / (m - n);
}

void PrintIteration(std::ofstream& outfile, int iteration, const std::vector<Point>& Points, const cv::Mat& X, const cv::Mat& V, const cv::Mat& L, const cv::Mat& So, int m, int n) {
    outfile << std::fixed << std::setprecision(3);

    outfile << "\n*************************************** Iteration " << iteration + 1 << " ***************************************\n\n";
    outfile << std::setw(10) << "So : " << sqrt(So.at<double>(0, 0)) << "\n\n";

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
}

void UpdateValues(std::vector<Point>& Points, const cv::Mat& X) {
    int n = Points.size() * 2;
    if (X.rows != n || X.cols != 1) {
        std::cerr << "Error: Size mismatch between X and Points vector." << std::endl;
        return;
    }

    for (int i = 0; i < n / 2; ++i) {
        Points[i].X += X.at<double>(2 * i, 0);
        Points[i].Y += X.at<double>(2 * i + 1, 0);
    }
}

void PrintResult(std::ofstream& outfile, int& m, const std::vector<Point>& Points, const std::vector<Point>& Control_Points, const std::vector<Line>& Distances, const std::vector<Angle>& Angles, const cv::Mat& X, const std::vector<Point>& Points_init, const std::vector<double>& predicteDistances, const std::vector<double>& predicteAngles, const cv::Mat& SigmaXX, const cv::Mat& SigmaLL, const cv::Mat& So, const double X2) {
    outfile << std::fixed << std::setprecision(3);

    outfile << "So : " << sqrt(So.at<double>(0, 0)) << "\n\n";

    int dof = m - Points.size() * 2;




    // 99% 신뢰 수준에서의 임계값 테이블
    std::map<int, std::pair<double, double>> chi_square_critical_values = {
        {1, {0.00004, 7.88}}, {2, {0.01, 10.60}}, {3, {0.07, 12.84}}, {4, {0.21, 14.86}},
        {5, {0.41, 16.75}}, {6, {0.68, 18.55}}, {7, {0.99, 20.28}}, {8, {1.34, 21.95}},
        {9, {1.73, 23.59}}, {10, {2.16, 25.19}}, {11, {2.60, 26.76}}, {12, {3.07, 28.30}},
        {13, {3.57, 29.82}}, {14, {4.07, 31.32}}, {15, {4.60, 32.80}}, {16, {5.14, 34.27}},
        {17, {5.70, 35.72}}, {18, {6.26, 37.16}}, {19, {6.84, 38.58}}, {20, {7.43, 40.00}},
        {21, {8.03, 41.40}}, {22, {8.64, 42.80}}, {23, {9.26, 44.18}}, {24, {9.89, 45.56}},
        {25, {10.52, 46.93}}, {26, {11.16, 48.29}}, {27, {11.81, 49.64}}, {28, {12.46, 50.99}},
        {29, {13.12, 52.34}}, {30, {13.79, 53.67}}, {40, {20.71, 63.69}}, {50, {27.99, 76.15}},
        {60, {35.53, 88.38}}, {70, {43.19, 100.42}}, {80, {51.17, 112.33}}, {90, {59.20, 124.12}},
        {100, {67.33, 135.81}}
    };

    // 입력받은 자유도에 대한 임계값을 찾기
    double critical_value_low = 0, critical_value_high = 0;

    outfile << "\n********************************** Chi-squared Test ***********************************\n";
    outfile << std::setw(5) <<" Degree of Freedom : " << dof << "\n";
	outfile << std::setw(5) << " X2 : " << X2 << "\n\n";

    if (chi_square_critical_values.find(dof) != chi_square_critical_values.end()) {
        critical_value_low = chi_square_critical_values[dof].first;
        critical_value_high = chi_square_critical_values[dof].second;
    }
    else {
        std::cerr << "No threshold was found for that degree of freedom.\n";
    }

	outfile << std::setw(5) << " [" << critical_value_low << " < " << X2 << " < " << critical_value_high << "]" << "\n";
    if (critical_value_low < X2 && X2 < critical_value_high) {
		outfile << "So^2 can be judged to be valid at a 99% confidence level.\n";
    }
    else {
		outfile << "It can be determined that So^2 is not valid at the 99% confidence level.\n";
    }


    outfile << "\n************************************* Point Data **************************************\n";
    outfile << std::setw(15) << "Point_ID" << "\t|"
        << std::setw(25) << "Point_Name" << "\t|"
        << std::setw(25) << "X_init" << "\t|"
        << std::setw(25) << "Y_init" << "\t|"
        << std::setw(25) << "X_final" << "\t|"
        << std::setw(25) << "Y_final" << "\t|"
        << std::setw(25) << "S_X" << "\t|"
        << std::setw(25) << "S_Y" << "\n\n";

    for (int i = 0; i < Points.size(); ++i) {
        outfile << std::setw(12) << Points[i].index << "\t"
            << std::setw(25) << Points[i].name << "\t"
            << std::setw(25) << Points_init[i].X << "\t"
            << std::setw(25) << Points_init[i].Y << "\t"
            << std::setw(25) << Points[i].X << "\t"
            << std::setw(25) << Points[i].Y << "\t"
            << std::setw(25) << sqrt(SigmaXX.at<double>(2 * i, 2 * i)) << "\t"
            << std::setw(25) << sqrt(SigmaXX.at<double>(2 * i + 1, 2 * i + 1)) << "\n\n";
    }

    outfile << "\n************************************* Line Data **************************************\n";
    outfile << std::setw(15) << "Line_ID" << "\t|"
        << std::setw(25) << "From" << "\t|"
        << std::setw(25) << "To" << "\t|"
        << std::setw(25) << "Distance_obs" << "\t|"
        << std::setw(25) << "Distance_cal" << "\t|"
        << std::setw(25) << "Distance_V" << "\t|"
		<< std::setw(25) << "S_Distance(m)" << "\n\n";

    for (int i = 0; i < Distances.size(); ++i) {
        outfile << std::setw(15) << Distances[i].index << "\t"
            << std::setw(25) << Distances[i].From << "\t"
            << std::setw(25) << Distances[i].To << "\t"
            << std::setw(25) << Distances[i].Distance << "\t"
            << std::setw(25) << predicteDistances[i] << "\t"
            << std::setw(25) << predicteDistances[i] - Distances[i].Distance << "\t"
            << std::setw(25) << sqrt(SigmaLL.at<double>(i, i)) << "\n\n";
    }

    outfile << "\n************************************* Angle Data **************************************\n";
    outfile << std::setw(15) << "Angle_ID" << "\t|"
        << std::setw(25) << "Backsight" << "\t|"
        << std::setw(25) << "Occupied" << "\t|"
        << std::setw(25) << "Foresight" << "\t|"
        << std::setw(25) << "Angle_obs" << "\t|"
        << std::setw(25) << "Angle_cal" << "\t|"
        << std::setw(25) << "Angle_V" << "\t|"
		<< std::setw(25) << "S_Angle(Deg)" << "\n\n";

    for (int i = 0; i < Angles.size(); ++i) {
        outfile << std::setw(15) << Angles[i].index << "\t"
            << std::setw(25) << Angles[i].Backsight << "\t"
            << std::setw(25) << Angles[i].Occupied << "\t"
            << std::setw(25) << Angles[i].Foresight << "\t"
            << std::setw(25) << Angles[i].Angle_D + Angles[i].Angle_M / 60.0 + Angles[i].Angle_S / 3600.0 << "\t"
            << std::setw(25) << predicteAngles[i] << "\t"
			<< std::setw(25) << predicteAngles[i] - (Angles[i].Angle_D + Angles[i].Angle_M / 60.0 + Angles[i].Angle_S / 3600.0) << "\t"
			<< std::setw(25) << sqrt(SigmaLL.at<double>(i + Distances.size(), i + Distances.size())*(180*M_PI)) << "\n\n";
    }





}
