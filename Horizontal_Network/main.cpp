#include "Horizontal_Network.h"
#include <QApplication>
#include <QMessageBox>
#include "function.h"


int main(int argc, char* argv[])
{


    QApplication a(argc, argv);
    Horizontal_Network w;
    w.show();

    // Horizontal_Network Ŭ������ �����͸� �������� ���� ���� ����
    QObject::connect(&w, &Horizontal_Network::startButtonClicked, [&w]() {
        std::string inputtitle, business_information, inputname, input_controldata, input_measuredata, input_angledata, input_linedata, output_data;


        /*********************************************************************
        1.  �Է� �ڷ� �б�(�۾����ü�)
        **********************************************************************/

        QString controlPath = w.getControlFilePath();
        QString measurePath = w.getMeasureFilePath();
        QString anglePath = w.getAngleFilePath();
        QString linePath = w.getLineFilePath();
        QString title = w.getTitle();
        QString business = w.getBusiness();
        QString name = w.getName();
        QString output = w.getOutput();

  //      inputtitle = title.toStdString();
  //      business_information = business.toStdString();
  //      inputname = name.toStdString();
  //      input_controldata = controlPath.toStdString();
		//input_measuredata = measurePath.toStdString();
		//input_angledata = anglePath.toStdString();
		//input_linedata = linePath.toStdString();
		//output_data = output.toStdString();

		inputtitle = "Horizontal Network Calculation";
		business_information = "Department of Civil En";
		inputname = "Kwon";
		input_controldata = "C:\\Horizontal_Network\\Input2\\C.csv";
		input_measuredata = "C:\\Horizontal_Network\\Input2\\M.csv";
		input_angledata = "C:\\Horizontal_Network\\Input2\\A.csv";
		input_linedata = "C:\\Horizontal_Network\\Input2\\D.csv";
		output_data = "C:\\Horizontal_Network\\Output\\Output.txt";



        std::vector<Point> points; // �Է� �ڷ�
		std::vector<Point> control_points; // ������
		std::vector<Point> measure_points; // ������
        std::vector<Point> points_init;
        std::vector<Line> lines;
        std::vector<Angle> angles;

        ReadPointFile(input_controldata, input_measuredata, points); // �Է� �ڷ� �б�
		ReadControlFile(input_controldata, control_points); // ������ �б�
		ReadMeasureFile(input_measuredata, measure_points); // ������ �б�
        ReadLineFile(input_linedata, lines);
        ReadAngleFile(input_angledata, angles);

        std::ofstream outfile(output_data);  // ��� �ڷ� ����

        PrintInformation(outfile, inputtitle, business_information, inputname, points, lines, angles); // �Է� �ڷ� ���


        /*********************************************************************
        2. ���� ��Ī ���İ� �� ���� �ο�
        **********************************************************************/
		SortFile(points, lines, angles); // ���� ��Ī ���İ� �� ���� �ο�
		PrintSortFile(outfile, points, lines, angles); // ���� ��Ī ���
		points_init = points;

        /*********************************************************************
        3. ����� ����
        **********************************************************************/
        outfile << "\n************************************* Start **************************************\n";

        int n = 2 * points.size(); // ������ ����
        int m = 2 * control_points.size() + lines.size() + angles.size(); // ���� ��
        cv::Mat W = cv::Mat::zeros(m, m, CV_64F); // W��� ����
		cv::Mat A = cv::Mat::zeros(m, n, CV_64F); // A��� ����
		cv::Mat L = cv::Mat::zeros(m, 1, CV_64F); // L��� ����
		cv::Mat X = cv::Mat::zeros(n, 1, CV_64F); // X��� ����
		cv::Mat V = cv::Mat::zeros(m, 1, CV_64F); // V��� ����
		cv::Mat So = cv::Mat::zeros(1, 1, CV_64F); // So��� ����
        std::vector<double> So_sqrt_vec;
        int iteration = 0;
        std::vector<double> predicteDistance;
		std::vector<double> predicteAngle;

        /*********************************************************************
        3.1 W ��� ����
        **********************************************************************/
		ComposeWMatrix(W, m, n, control_points, lines, angles); // W��� ����

        /*********************************************************************
        3.2 �ݺ� ����
        **********************************************************************/
        for (iteration; iteration < 100; ++iteration) {
			ComposeAMatrix(A, m, n, points, lines, angles, control_points); // A��� ����
			ComposeLMatrix(L, m, points, lines, angles, control_points, predicteDistance, predicteAngle); // L��� ����
			CalculateXMatrix(X, A, W, L, m, n); // X��� ���
			CalculateVMatrix(V, A, X, L, m, n); // V��� ���
			CalculateSoMatrix(So, V, W, X, m, n); // So��� ���
			So_sqrt_vec.push_back(std::sqrt(So.at<double>(0, 0))); // So�� �������� ���Ϳ� ����

			PrintIteration(outfile, iteration, points, X, V, L , So, m, n); // �ݺ� ���

            double max_diff = 0.0;
            for (int i = 0; i < X.rows; ++i)
            {
				max_diff = std::max(max_diff, std::abs(X.at<double>(i, 0)));
            }

            if (iteration == 100)
            {
                outfile << "\n************************************* End **************************************\n";
                outfile << "Iteration : " << iteration+1 << "  Adjustment termination : Maximum number of repeated calculations (abnormal termination)\n";
                UpdateValues(points, X, n); // ���� ��ǥ ������Ʈ


                break;
            }

            else if (max_diff < 0.000001) {
                outfile << "\n************************************* End **************************************\n";
                outfile << "Iteration : " << iteration + 1 << "  Adjustment termination: When the largest adjustment value is reduced to a certain extent, it is terminated (normal termination)n";
                UpdateValues(points, X, n); // ���� ��ǥ ������Ʈ

                break;
            }

			else if (So_sqrt_vec.size() > 1 && std::abs(So_sqrt_vec[So_sqrt_vec.size() - 1] - So_sqrt_vec[So_sqrt_vec.size() - 2]) < 0.000001)
            {
				outfile << "\n************************************* End **************************************\n";
				outfile << "Iteration : " << iteration + 1 << "  Adjustment termination: When the adjustment value is reduced to a certain extent, it is terminated (normal termination)\n";
                UpdateValues(points, X, n); // ���� ��ǥ ������Ʈ

				break;
			}

			else if (iteration >= 2 && std::abs(So_sqrt_vec[iteration] - So_sqrt_vec[iteration - 1]) > 0)
            {
                if (So_sqrt_vec[iteration - 1] > So_sqrt_vec[iteration - 2]) {
					outfile << "\n************************************* End **************************************\n";
					outfile << "Iteration : " << iteration + 1 << "  Adjustment termination: When the adjustment value is increased, it is terminated (abnormal termination)\n";
                    UpdateValues(points, X, n); // ���� ��ǥ ������Ʈ

					break;
                }
			}

			UpdateValues(points, X, n); // ���� ��ǥ ������Ʈ
        }

		PrintResult(outfile, points, control_points, lines, angles, X, points_init, n, predicteDistance, predicteAngle); // ��� ���






        outfile.close();


        });



    return a.exec();
}
