#include "Horizontal_Network.h"
#include <QApplication>
#include <QMessageBox>
#include "function.h"


int main(int argc, char* argv[])
{


    QApplication a(argc, argv);
    Horizontal_Network w;
    w.show();

    // Horizontal_Network 클래스의 데이터를 가져오기 위한 슬롯 연결
    QObject::connect(&w, &Horizontal_Network::startButtonClicked, [&w]() {
        std::string inputtitle, business_information, inputname, input_controldata, input_measuredata, input_angledata, input_linedata, output_data;

        /*********************************************************************
        1.  입력 자료 읽기(작업지시서)
        **********************************************************************/

        QString controlPath = w.getControlFilePath();
        QString measurePath = w.getMeasureFilePath();
        QString anglePath = w.getAngleFilePath();
        QString linePath = w.getLineFilePath();
        QString title = w.getTitle();
        QString business = w.getBusiness();
        QString name = w.getName();
        QString output = w.getOutput();

        inputtitle = title.toStdString();
        business_information = business.toStdString();
        inputname = name.toStdString();
        input_controldata = controlPath.toStdString();
		input_measuredata = measurePath.toStdString();
		input_angledata = anglePath.toStdString();
		input_linedata = linePath.toStdString();
		output_data = output.toStdString();


        std::vector<Point> points; // 입력 자료
		std::vector<Point> control_points; // 제어점
		std::vector<Point> measure_points; // 측정점
        std::vector<Point> points_init;
        std::vector<Line> lines;
        std::vector<Angle> angles;

        ReadPointFile(input_controldata, input_measuredata, points); // 입력 자료 읽기
		ReadControlFile(input_controldata, control_points); // 제어점 읽기
		ReadMeasureFile(input_measuredata, measure_points); // 측정점 읽기
        ReadLineFile(input_linedata, lines); // 거리 자료 읽기
        ReadAngleFile(input_angledata, angles); // 각 자료 읽기

        std::ofstream outfile(output_data);  // 출력 자료 생성

        PrintInformation(outfile, inputtitle, business_information, inputname, points, lines, angles); // 입력 자료 출력


        /*********************************************************************
        2. 측점 명칭 정렬과 논리 순서 부여
        **********************************************************************/
		SortFile(points, lines, angles); // 측점 명칭 정렬과 논리 순서 부여
		PrintSortFile(outfile, points, lines, angles); // 측점 명칭 출력
		points_init = points;

        /*********************************************************************
        3. 수평망 조정
        **********************************************************************/
        outfile << "\n****************************************** Start ******************************************\n";

        int n = 2 * points.size(); // 측점의 개수
        int m = 2 * control_points.size() + lines.size() + angles.size(); // 관측 수
        cv::Mat W = cv::Mat::zeros(m, m, CV_64F); // W행렬 생성
		cv::Mat A = cv::Mat::zeros(m, n, CV_64F); // A행렬 생성
		cv::Mat L = cv::Mat::zeros(m, 1, CV_64F); // L행렬 생성
		cv::Mat X = cv::Mat::zeros(n, 1, CV_64F); // X행렬 생성
		cv::Mat V = cv::Mat::zeros(m, 1, CV_64F); // V행렬 생성
		cv::Mat So = cv::Mat::zeros(1, 1, CV_64F); // So행렬 생성
		cv::Mat SigmaXX = cv::Mat::zeros(n, n, CV_64F); // SigmaXX행렬 생성
		cv::Mat SigmaLL = cv::Mat::zeros(m, m, CV_64F); // SigmaLL행렬 생성
        std::vector<double> predicteDistance;
        std::vector<double> predicteAngle;
        std::vector<double> So_sqrt_vec;
        int iteration = 0;


        /*********************************************************************
        3.1 W 행렬 조성
        **********************************************************************/
		ComposeWMatrix(W, control_points, lines, angles); // W행렬 조성


        /*********************************************************************
        3.2 반복 시작
        **********************************************************************/
        for (iteration; iteration < 100; ++iteration) {
            predicteDistance.clear();
            predicteAngle.clear();

            ComposeAMatrix(A, points, lines, angles, control_points); // A행렬 조성
            ComposeLMatrix(L, points, lines, angles, control_points, predicteDistance, predicteAngle); // L행렬 조성
            CalculateXMatrix(X, A, W, L); // X행렬 계산
            CalculateVMatrix(V, A, X, L); // V행렬 계산
            CalculateSoMatrix(So, V, W, m, n); // So행렬 계산
            So_sqrt_vec.push_back(std::sqrt(So.at<double>(0, 0))); // So의 제곱근을 벡터에 저장
            SigmaXX = So.at<double>(0, 0) * (A.t() * W * A).inv(); // SigmaXX 계산
			SigmaLL = A * SigmaXX * A.t(); // SigmaLL 계산
            PrintIteration(outfile, iteration, points, X, V, L, So, m, n); // 반복 출력

            double max_diff = 0.0;
            for (int i = 0; i < X.rows; ++i) {
                max_diff = std::max(max_diff, std::abs(X.at<double>(i, 0)));
            }

            if (iteration == 99) {
                outfile << "\n***************************************** End *****************************************\n";
                outfile << "Iteration : " << iteration + 1 << "  \nAdjustment termination : Maximum number of repeated calculations (abnormal termination)\n";
                UpdateValues(points, X); // 측점 좌표 업데이트
                break;
            }

            else if (max_diff < 0.000001) {
                outfile << "\n***************************************** End *****************************************\n";
                outfile << "Iteration : " << iteration + 1 << "  \nAdjustment termination: When the largest adjustment value is reduced to a certain extent, it is terminated (normal termination)\n";
                UpdateValues(points, X); // 측점 좌표 업데이트
                break;
            }

            else if (So_sqrt_vec.size() > 1 && std::abs(So_sqrt_vec[So_sqrt_vec.size() - 1] - So_sqrt_vec[So_sqrt_vec.size() - 2]) < 0.000001) {
                outfile << "\n***************************************** End *****************************************\n";
                outfile << "Iteration : " << iteration + 1 << "  \nAdjustment termination: When the adjustment value is reduced to a certain extent, it is terminated (normal termination)\n";
                UpdateValues(points, X); // 측점 좌표 업데이트
                break;
            }

            else if (iteration >= 2 && std::abs(So_sqrt_vec[iteration] - So_sqrt_vec[iteration - 1]) > 0) {
                if (So_sqrt_vec[iteration - 1] > So_sqrt_vec[iteration - 2]) {
                    outfile << "\n***************************************** End *****************************************\n";
                    outfile << "Iteration : " << iteration + 1 << "  \nAdjustment termination: When the adjustment value is increased, it is terminated (abnormal termination)\n";
                    UpdateValues(points, X); // 측점 좌표 업데이트
                    break;
                }
            }

            UpdateValues(points, X); // 측점 좌표 업데이트
        }

		double X_test = (m - n) * sqrt(So.at<double>(0, 0)) * sqrt(So.at<double>(0, 0)) / (So_sqrt_vec[iteration - 1] * So_sqrt_vec[iteration - 1]);
		PrintResult(outfile, m, points, control_points, lines, angles, X, points_init, predicteDistance, predicteAngle, SigmaXX, SigmaLL, So, X_test); // 결과 출력

		outfile.close(); // 출력 자료 닫기
        });



    return a.exec();
}
