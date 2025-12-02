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
    vector<double> decodef(const string& str, int n) {
        vector<double> result(n, 0.0);
        istringstream iss(str);
        string token;
        int i = 0;
        
        while (iss >> token && i < n) {
            result[i++] = stod(token);
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
                
                // 设计矩阵
                H(i, 0) = -dx/range;
                H(i, 1) = -dy/range;
                H(i, 2) = -dz/range;
                H(i, 3) = 1.0;  // 钟差参数
                
                // 残差 = 伪距观测值 - 几何距离 - 钟差
                res(i) = obs[i].pseudorange - range - x[3];
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
            
            if (max_dx < 1e-4) {
                cout << "Converged after " << iter + 1 << " iterations" << endl;
                break;
            }
        }
        
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
        fout << "GPSS Week    Seconds of Week          X(m)              Y(m)              Z(m)           Clock(m)" << endl;
        
        string line;
        int epochnum, gpsw, obsnum = 0;
        double gpss;
        vector<SatelliteData> currentEpoch;
        
        while (getline(fin, line)) {
            if (line.length() < 20) continue;
            
            if (line[0] == '#') {
                // 如果是新历元的开始，先处理上一个历元的数据
                if (!currentEpoch.empty()) {
                    vector<double> x;
                    if (pntpos(currentEpoch, x)) {
                        fout << setw(4) << gpsw << " " 
                             << fixed << setprecision(3) << setw(14) << gpss << " "
                             << setw(14) << x[0] << " " << setw(14) << x[1] << " " 
                             << setw(14) << x[2] << " " << setw(14) << x[3] << endl;
                    } else {
                        fout << setw(4) << gpsw << " " 
                             << fixed << setprecision(3) << setw(14) << gpss << " "
                             << setw(14) << 0.0 << " " << setw(14) << 0.0 << " " 
                             << setw(14) << 0.0 << " " << setw(14) << 0.0 << endl;
                    }
                    currentEpoch.clear();
                }
                
                // 解析新历元头信息
                istringstream iss(line.substr(2));
                iss >> epochnum >> gpsw >> gpss >> obsnum;
                cout << "Epoch: " << epochnum << ", GPS Week: " << gpsw 
                     << ", GPS Seconds: " << gpss << ", Obs: " << obsnum << endl;
                
            } else if (obsnum > 0 && currentEpoch.size() < obsnum && currentEpoch.size() < MAXOBS) {
                // 解析卫星数据
                string satid = line.substr(0, 3);
                vector<double> values = decodef(line.substr(3), 4);
                
                if (values.size() == 4) {
                    currentEpoch.emplace_back(satid, values[0], values[1], values[2], values[3]);
                }
            }
        }
        
        // 处理最后一个历元
        if (!currentEpoch.empty()) {
            vector<double> x;
            if (pntpos(currentEpoch, x)) {
                fout << setw(4) << gpsw << " " 
                     << fixed << setprecision(3) << setw(14) << gpss << " "
                     << setw(14) << x[0] << " " << setw(14) << x[1] << " " 
                     << setw(14) << x[2] << " " << setw(14) << x[3] << endl;
            } else {
                fout << setw(4) << gpsw << " " 
                     << fixed << setprecision(3) << setw(14) << gpss << " "
                     << setw(14) << 0.0 << " " << setw(14) << 0.0 << " " 
                     << setw(14) << 0.0 << " " << setw(14) << 0.0 << endl;
            }
        }
        
        fin.close();
        fout.close();
        cout << "Processing completed. Results saved to " << outputFile << endl;
    }
};

int main() {
    PointPositioning processor;
    
    string inputFile = "E:/STUDY/Sophomore1/最优估计/第二次上机实习/work2/2059329.25ObsCorr";
    string outputFile = "outpos.txt";
    
    processor.processFile(inputFile, outputFile);
    
    return 0;
}