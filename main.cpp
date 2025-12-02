#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const int MAXOBS = 50;

class SatelliteData {
public:
    string satid;
    double X, Y, Z, pseudorange;
    
    SatelliteData(const string& id, double x, double y, double z, double pr) 
        : satid(id), X(x), Y(y), Z(z), pseudorange(pr) {}
};

class PointPositioning {
private:
    vector<double> decodeLine(const string& line) {
        vector<double> result;
        istringstream iss(line);
        string token;
        
        // 跳过卫星PRN号（前3个字符）
        string satid;
        iss >> satid;
        
        // 读取后面的4个数值：X, Y, Z, pseudorange
        for (int i = 0; i < 4; i++) {
            if (iss >> token) {
                result.push_back(stod(token));
            }
        }
        
        return result;
    }
    
public:
    bool pntpos(const vector<SatelliteData>& obs, vector<double>& x) {
        if (obs.size() < 4) {
            cout << "Not enough satellites (" << obs.size() << " < 4)" << endl;
            return false;
        }
        
        int n = obs.size();
        
        // 初始位置设为地球中心
        x = vector<double>(4, 0.0); // x, y, z, clock
        
        cout << "Processing " << n << " satellites: ";
        for (const auto& sat : obs) {
            cout << sat.satid << " ";
        }
        cout << endl;
        
        int iter;
        for (iter = 0; iter < 10; iter++) {
            // 构建设计矩阵和残差向量
            MatrixXd H(n, 4);
            VectorXd res(n);
            
            for (int i = 0; i < n; i++) {
                double dx = obs[i].X - x[0];
                double dy = obs[i].Y - x[1];
                double dz = obs[i].Z - x[2];
                double range = sqrt(dx*dx + dy*dy + dz*dz);
                
                if (range < 1e-12) {
                    cout << "Zero range for satellite " << obs[i].satid << endl;
                    return false;
                }
                
                // 设计矩阵
                H(i, 0) = -dx/range;
                H(i, 1) = -dy/range;
                H(i, 2) = -dz/range;
                H(i, 3) = 1.0;  // 钟差参数
                
                // 残差 = 伪距观测值 - 几何距离 - 钟差
                res(i) = obs[i].pseudorange - range - x[3];
                
                if (iter == 0) {
                    cout << "Sat " << obs[i].satid << ": range=" << range 
                         << ", pseudorange=" << obs[i].pseudorange 
                         << ", residual=" << res(i) << endl;
                }
            }
            
            // 权重矩阵（单位阵）
            MatrixXd W = MatrixXd::Identity(n, n);
            
            // 法方程：H^T * W * H * dx = H^T * W * res
            MatrixXd HTH = H.transpose() * W * H;
            VectorXd HTR = H.transpose() * W * res;
            
            // 检查矩阵是否可逆
            FullPivLU<MatrixXd> lu(HTH);
            if (!lu.isInvertible()) {
                cout << "Singular matrix at iteration " << iter << endl;
                return false;
            }
            
            // 解法方程
            VectorXd dx_vec = HTH.ldlt().solve(HTR);
            
            // 更新参数
            x[0] += dx_vec(0);
            x[1] += dx_vec(1);
            x[2] += dx_vec(2);
            x[3] += dx_vec(3);
            
            // 检查收敛
            double max_dx = 0.0;
            for (int i = 0; i < 3; i++) {
                if (abs(dx_vec(i)) > max_dx) {
                    max_dx = abs(dx_vec(i));
                }
            }
            
            cout << "Iteration " << iter + 1 << ": dx=(" << dx_vec(0) << ", " 
                 << dx_vec(1) << ", " << dx_vec(2) << "), max_dx=" << max_dx << endl;
            
            if (max_dx < 1e-4) {
                cout << "Converged after " << iter + 1 << " iterations" << endl;
                break;
            }
        }
        
        cout << "Final position: (" << x[0] << ", " << x[1] << ", " << x[2] 
             << "), clock=" << x[3] << endl;
        
        return true;
    }
    
    void processFile(const string& inputFile, const string& outputFile) {
        ifstream fin(inputFile);
        ofstream fout(outputFile);
        
        if (!fin.is_open()) {
            cerr << "Error opening input file: " << inputFile << endl;
            return;
        }
        
        if (!fout.is_open()) {
            cerr << "Error opening output file: " << outputFile << endl;
            fin.close();
            return;
        }
        
        // 写入表头
        fout << "GPSWeek  SecondsOfWeek          X(m)              Y(m)              Z(m)           Clock(m)  NumSats" << endl;
        
        string line;
        int epochnum, gpsw, obsnum = 0;
        double gpss;
        vector<SatelliteData> currentEpoch;
        
        int epochCount = 0;
        
        while (getline(fin, line)) {
            if (line.length() < 10) continue;
            
            if (line[0] == '#') {
                // 如果是新历元的开始，先处理上一个历元的数据
                if (!currentEpoch.empty()) {
                    epochCount++;
                    cout << "\n=== Processing Epoch " << epochCount << " ===" << endl;
                    
                    vector<double> x;
                    if (pntpos(currentEpoch, x)) {
                        fout << setw(4) << gpsw << " " 
                             << fixed << setprecision(3) << setw(14) << gpss << " "
                             << setw(14) << x[0] << " " << setw(14) << x[1] << " " 
                             << setw(14) << x[2] << " " << setw(14) << x[3] << " "
                             << setw(6) << currentEpoch.size() << endl;
                    } else {
                        fout << setw(4) << gpsw << " " 
                             << fixed << setprecision(3) << setw(14) << gpss << " "
                             << setw(14) << 0.0 << " " << setw(14) << 0.0 << " " 
                             << setw(14) << 0.0 << " " << setw(14) << 0.0 << " "
                             << setw(6) << currentEpoch.size() << endl;
                    }
                    currentEpoch.clear();
                }
                
                // 解析新历元头信息
                istringstream iss(line.substr(2));
                iss >> epochnum >> gpsw >> gpss >> obsnum;
                cout << "Epoch: " << epochnum << ", GPS Week: " << gpsw 
                     << ", GPS Seconds: " << gpss << ", Expected Obs: " << obsnum << endl;
                
            } else if (obsnum > 0 && currentEpoch.size() < obsnum && currentEpoch.size() < MAXOBS) {
                // 解析卫星数据
                string satid = line.substr(0, 3);
                vector<double> values = decodeLine(line);
                
                if (values.size() == 4) {
                    currentEpoch.emplace_back(satid, values[0], values[1], values[2], values[3]);
                    cout << "Added satellite: " << satid << " with pseudorange: " << values[3] << endl;
                } else {
                    cout << "Warning: Invalid data line for satellite " << satid 
                         << ", expected 4 values, got " << values.size() << endl;
                }
            }
        }
        
        // 处理最后一个历元
        if (!currentEpoch.empty()) {
            epochCount++;
            cout << "\n=== Processing Epoch " << epochCount << " ===" << endl;
            
            vector<double> x;
            if (pntpos(currentEpoch, x)) {
                fout << setw(4) << gpsw << " " 
                     << fixed << setprecision(3) << setw(14) << gpss << " "
                     << setw(14) << x[0] << " " << setw(14) << x[1] << " " 
                     << setw(14) << x[2] << " " << setw(14) << x[3] << " "
                     << setw(6) << currentEpoch.size() << endl;
            } else {
                fout << setw(4) << gpsw << " " 
                     << fixed << setprecision(3) << setw(14) << gpss << " "
                     << setw(14) << 0.0 << " " << setw(14) << 0.0 << " " 
                     << setw(14) << 0.0 << " " << setw(14) << 0.0 << " "
                     << setw(6) << currentEpoch.size() << endl;
            }
        }
        
        fin.close();
        fout.close();
        cout << "\nProcessing completed. Results saved to " << outputFile << endl;
    }
};

int main() {
    PointPositioning processor;
    
    string inputFile = "E:/STUDY/Sophomore1/最优估计/第二次上机实习/work2/2059329.25ObsCorr";  // 修改为您的实际文件路径
    string outputFile = "outpos.txt";
    
    processor.processFile(inputFile, outputFile);
    
    return 0;
}