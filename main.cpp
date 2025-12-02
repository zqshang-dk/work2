#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// 常数定义
const double PI = 3.14159265358979323846;
const double C = 299792458.0;  // 光速，m/s
const double OMEGA_E = 7.2921151467e-5;  // 地球自转角速度，rad/s
const double RE_WGS84 = 6378137.0;  // WGS84椭球长半轴
const double FE_WGS84 = 1.0 / 298.257223563;  // WGS84扁率

// 结构体定义
struct Observation {
    string satID;      // 卫星ID
    double X, Y, Z;    // 卫星坐标，m
    double pseudorange; // 伪距观测值，m
    double elevation;   // 高度角，度（如果有）
};

struct Epoch {
    int epochNum;      // 历元编号
    int gpsWeek;       // GPS周
    double gpsSeconds; // GPS秒
    int obsNum;        // 观测卫星数
    vector<Observation> obs; // 观测数据
};

struct PositionResult {
    int gpsWeek;
    double gpsSeconds;
    double X, Y, Z, T; // 接收机坐标和钟差
    double sigma2;     // 验后单位权方差
    MatrixXd Qxx;      // 协方差矩阵
    Vector3d ENU;      // ENU坐标
};

// 函数声明
vector<Epoch> readObservationFile(const string& filename);
PositionResult pointPositioning(const vector<Observation>& obs, const Vector4d& approx, bool firstEpoch);
void writeResults(const string& filename, const vector<PositionResult>& results);
void convertXYZtoENU(double X, double Y, double Z, 
                     double refX, double refY, double refZ,
                     double& E, double& N, double& U);
void calculateStatistics(const vector<PositionResult>& results, 
                        double refX, double refY, double refZ);
Vector2d calcLL(double X, double Y, double Z);

// 主函数
int main() {
    // 1. 读取观测数据文件
    string inputFile = "E:/STUDY/Sophomore1/最优估计/第二次上机实习/work2/CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt";
    string outputFile = "position_results.txt";
    
    cout << "正在读取观测数据文件: " << inputFile << endl;
    vector<Epoch> epochs = readObservationFile(inputFile);
    cout << "成功读取 " << epochs.size() << " 个历元的数据" << endl;
    
    if (epochs.empty()) {
        cerr << "错误：没有读取到有效数据！" << endl;
        return -1;
    }
    
    // 2. 逐历元进行定位解算
    vector<PositionResult> results;
    Vector4d approx = Vector4d::Zero();  // 初始近似值 [X, Y, Z, T]
    
    for (size_t i = 0; i < epochs.size(); i++) {
        const Epoch& epoch = epochs[i];
        
        // 检查是否有足够卫星
        if (epoch.obsNum < 4) {
            cerr << "历元 " << epoch.epochNum << " 卫星数不足4颗，跳过" << endl;
            continue;
        }
        
        cout << "处理历元 " << epoch.epochNum 
             << " (" << epoch.gpsWeek << " " << epoch.gpsSeconds 
             << "), 卫星数: " << epoch.obsNum << endl;
        
        // 进行单历元定位解算
        PositionResult result = pointPositioning(epoch.obs, approx, (i == 0));
        
        // 更新近似值供下一个历元使用
        approx(0) = result.X;
        approx(1) = result.Y;
        approx(2) = result.Z;
        approx(3) = result.T;
        
        // 设置历元信息
        result.gpsWeek = epoch.gpsWeek;
        result.gpsSeconds = epoch.gpsSeconds;
        
        // 存储结果
        results.push_back(result);
        
        // 显示进度
        if ((i + 1) % 100 == 0 || i == epochs.size() - 1) {
            cout << "已处理 " << i + 1 << "/" << epochs.size() << " 个历元" << endl;
        }
    }
    
    cout << "定位解算完成，共处理 " << results.size() << " 个历元" << endl;
    
    // 3. 写入结果文件
    writeResults(outputFile, results);
    cout << "结果已写入: " << outputFile << endl;
    
    // 4. 计算统计信息
    // 参考坐标（任务书中给出）
    double refX = -2148744.3960;
    double refY = 4426641.4191;
    double refZ = 4044655.5364;
    
    calculateStatistics(results, refX, refY, refZ);
    
    // 5. 生成绘图数据文件（用于MATLAB/Python绘图）
    ofstream enuFile("enu_results.txt");
    if (enuFile.is_open()) {
        enuFile << "Epoch GPSWeek GPSSeconds E N U" << endl;
        for (size_t i = 0; i < results.size(); i++) {
            enuFile << i + 1 << " " 
                   << results[i].gpsWeek << " " 
                   << fixed << setprecision(3) << results[i].gpsSeconds << " "
                   << scientific << setprecision(6)
                   << results[i].ENU(0) << " "
                   << results[i].ENU(1) << " "
                   << results[i].ENU(2) << endl;
        }
        enuFile.close();
        cout << "ENU坐标已写入: enu_results.txt" << endl;
    }
    
    return 0;
}

// 读取观测数据文件
vector<Epoch> readObservationFile(const string& filename) {
    vector<Epoch> epochs;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "错误：无法打开文件 " << filename << endl;
        return epochs;
    }
    
    string line;
    Epoch currentEpoch;
    bool readingObservations = false;
    int obsCount = 0;
    
    while (getline(file, line)) {
        // 跳过空行
        if (line.empty()) continue;
        
        // 检查是否为历元头
        if (line[0] == '#') {
            // 如果有上一个历元的数据，先保存
            if (readingObservations && currentEpoch.obsNum > 0) {
                epochs.push_back(currentEpoch);
            }
            
            // 解析新的历元头
            istringstream iss(line.substr(2)); // 跳过"# "
            iss >> currentEpoch.epochNum >> currentEpoch.gpsWeek 
                >> currentEpoch.gpsSeconds >> currentEpoch.obsNum;
            
            currentEpoch.obs.clear();
            readingObservations = true;
            obsCount = 0;
        } 
        else if (readingObservations && obsCount < currentEpoch.obsNum) {
            // 解析观测数据行
            Observation obs;
            istringstream iss(line);
            iss >> obs.satID >> obs.X >> obs.Y >> obs.Z 
                >> obs.pseudorange >> obs.elevation;
            
            currentEpoch.obs.push_back(obs);
            obsCount++;
            
            // 如果已读完所有观测，保存历元
            if (obsCount == currentEpoch.obsNum) {
                epochs.push_back(currentEpoch);
                readingObservations = false;
            }
        }
    }
    
    // 处理最后一个历元
    if (readingObservations && currentEpoch.obsNum > 0) {
        epochs.push_back(currentEpoch);
    }
    
    file.close();
    return epochs;
}

// 单历元最小二乘定位解算
PositionResult pointPositioning(const vector<Observation>& obs, 
                                const Vector4d& approx, bool firstEpoch) {
    PositionResult result;
    int nObs = obs.size();
    
    // 迭代解算
    Vector4d X = approx;
    Vector4d dX;
    Matrix4d Qxx;
    double sigma2 = 0;
    
    int maxIter = 10;
    double tolerance = 1e-4;
    bool converged = false;
    
    for (int iter = 0; iter < maxIter; iter++) {
        // 构建设计矩阵G和观测值向量L
        MatrixXd G(nObs, 4);
        VectorXd L(nObs);
        VectorXd W(nObs);  // 权矩阵对角线元素
        
        for (int i = 0; i < nObs; i++) {
            const Observation& o = obs[i];
            
            // 计算几何距离
            double dx = o.X - X(0);
            double dy = o.Y - X(1);
            double dz = o.Z - X(2);
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            // 设计矩阵行
            G(i, 0) = -dx / r;
            G(i, 1) = -dy / r;
            G(i, 2) = -dz / r;
            G(i, 3) = 1.0;  // 钟差参数
            
            // 观测值向量（伪距残差）
            L(i) = o.pseudorange - (r + X(3) * C);
            
            // 权值（基于高度角，简单模型）
            double elev_rad = o.elevation * PI / 180.0;
            W(i) = sin(elev_rad) * sin(elev_rad);  // 高度角越大，权值越大
        }
        
        // 最小二乘解算：X = (G^T * W * G)^-1 * G^T * W * L
        MatrixXd W_mat = W.asDiagonal();
        Matrix4d GTWG = G.transpose() * W_mat * G;
        
        // 使用Eigen的矩阵求逆
        if (GTWG.determinant() == 0) {
            cerr << "错误：法方程矩阵奇异！" << endl;
            break;
        }
        
        dX = GTWG.inverse() * G.transpose() * W_mat * L;
        
        // 更新参数估计
        X += dX;
        
        // 检查收敛性
        if (dX.norm() < tolerance) {
            converged = true;
            
            // 计算验后单位权方差
            VectorXd V = L - G * dX;
            double VWV = V.transpose() * W_mat * V;
            sigma2 = VWV / (nObs - 4);
            
            // 计算协方差矩阵
            Qxx = sigma2 * GTWG.inverse();
            
            break;
        }
    }
    
    if (!converged) {
        cerr << "警告：迭代未收敛" << endl;
    }
    
    // 存储结果
    result.X = X(0);
    result.Y = X(1);
    result.Z = X(2);
    result.T = X(3);
    result.sigma2 = sigma2;
    result.Qxx = Qxx;
    
    return result;
}

// 写入结果文件
void writeResults(const string& filename, const vector<PositionResult>& results) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "错误：无法创建输出文件 " << filename << endl;
        return;
    }
    
    // 写入表头
    file << "GPSWeek  GPSSec             X(m)             Y(m)             Z(m)           Clock(m)     Sigma2" << endl;
    file << "=================================================================================================" << endl;
    
    // 写入数据
    for (const auto& result : results) {
        file << fixed << setw(7) << result.gpsWeek << " "
             << setw(12) << setprecision(3) << result.gpsSeconds << " "
             << scientific << setprecision(6)
             << setw(16) << result.X << " "
             << setw(16) << result.Y << " "
             << setw(16) << result.Z << " "
             << setw(16) << result.T * C << " "  // 钟差转换为米
             << setw(12) << result.sigma2 << endl;
    }
    
    file.close();
}

// 计算经纬度
Vector2d calcLL(double X, double Y, double Z) {
    double a = RE_WGS84;
    double f = FE_WGS84;
    double e2 = 2*f - f*f;
    
    double p = sqrt(X*X + Y*Y);
    double theta = atan2(Z * a, p * (1 - e2));
    
    double lat = atan2(Z + e2 * (1 - e2) * a * pow(sin(theta), 3) / (1 - e2),
                       p - e2 * a * pow(cos(theta), 3));
    double lon = atan2(Y, X);
    
    return Vector2d(lat, lon);
}

// XYZ转ENU坐标
void convertXYZtoENU(double X, double Y, double Z,
                     double refX, double refY, double refZ,
                     double& E, double& N, double& U) {
    
    // 计算参考点的经纬度
    Vector2d refLL = calcLL(refX, refY, refZ);
    double refLat = refLL(0);
    double refLon = refLL(1);
    
    // 计算旋转矩阵
    double sinLat = sin(refLat);
    double cosLat = cos(refLat);
    double sinLon = sin(refLon);
    double cosLon = cos(refLon);
    
    // 计算ENU
    double dX = X - refX;
    double dY = Y - refY;
    double dZ = Z - refZ;
    
    E = -sinLon * dX + cosLon * dY;
    N = -sinLat * cosLon * dX - sinLat * sinLon * dY + cosLat * dZ;
    U = cosLat * cosLon * dX + cosLat * sinLon * dY + sinLat * dZ;
}

// 计算统计信息
void calculateStatistics(const vector<PositionResult>& results,
                        double refX, double refY, double refZ) {
    
    if (results.empty()) return;
    
    int n = results.size();
    double sumE = 0, sumN = 0, sumU = 0;
    double sumE2 = 0, sumN2 = 0, sumU2 = 0;
    
    cout << "\n========== 统计结果 ==========" << endl;
    
    // 创建非const的results引用以便修改ENU
    vector<PositionResult>& mutableResults = const_cast<vector<PositionResult>&>(results);
    
    // 计算ENU坐标和统计信息
    for (int i = 0; i < n; i++) {
        double E, N, U;
        convertXYZtoENU(results[i].X, results[i].Y, results[i].Z,
                        refX, refY, refZ, E, N, U);
        
        // 存储ENU坐标
        // results[i].ENU << E, N, U;
        
        // 正确的方式：逐个分量赋值给ENU向量
        mutableResults[i].ENU(0) = E;
        mutableResults[i].ENU(1) = N;
        mutableResults[i].ENU(2) = U;


        sumE += E;
        sumN += N;
        sumU += U;
        sumE2 += E * E;
        sumN2 += N * N;
        sumU2 += U * U;
    }
    
    double meanE = sumE / n;
    double meanN = sumN / n;
    double meanU = sumU / n;
    
    double rmsE = sqrt(sumE2 / n);
    double rmsN = sqrt(sumN2 / n);
    double rmsU = sqrt(sumU2 / n);
    
    // 输出统计信息
    cout << fixed << setprecision(4);
    cout << "参考坐标: X=" << refX << " Y=" << refY << " Z=" << refZ << endl;
    cout << "历元总数: " << n << endl;
    cout << "\nENU坐标统计:" << endl;
    cout << "方向   均值(m)     RMS(m)" << endl;
    cout << "-------------------------" << endl;
    cout << "东(E)  " << setw(9) << meanE << "  " << setw(9) << rmsE << endl;
    cout << "北(N)  " << setw(9) << meanN << "  " << setw(9) << rmsN << endl;
    cout << "天顶(U)" << setw(9) << meanU << "  " << setw(9) << rmsU << endl;
    
    // 计算XYZ坐标统计
    double meanX = 0, meanY = 0, meanZ = 0;
    double stdX = 0, stdY = 0, stdZ = 0;
    
    for (int i = 0; i < n; i++) {
        meanX += results[i].X;
        meanY += results[i].Y;
        meanZ += results[i].Z;
    }
    meanX /= n; meanY /= n; meanZ /= n;
    
    for (int i = 0; i < n; i++) {
        stdX += pow(results[i].X - meanX, 2);
        stdY += pow(results[i].Y - meanY, 2);
        stdZ += pow(results[i].Z - meanZ, 2);
    }
    stdX = sqrt(stdX / (n - 1));
    stdY = sqrt(stdY / (n - 1));
    stdZ = sqrt(stdZ / (n - 1));
    
    cout << "\nXYZ坐标统计:" << endl;
    cout << "均值: X=" << scientific << setprecision(6) << meanX 
         << " Y=" << meanY << " Z=" << meanZ << endl;
    cout << fixed << setprecision(4);
    cout << "标准差: StdX=" << stdX << " StdY=" << stdY << " StdZ=" << stdZ << endl;
    
    // 输出精度信息
    cout << "\n精度评定 (1-sigma):" << endl;
    if (!results.empty()) {
        const PositionResult& lastResult = results.back();
        cout << "X方向精度: " << sqrt(lastResult.Qxx(0,0)) << " m" << endl;
        cout << "Y方向精度: " << sqrt(lastResult.Qxx(1,1)) << " m" << endl;
        cout << "Z方向精度: " << sqrt(lastResult.Qxx(2,2)) << " m" << endl;
        cout << "钟差精度: " << sqrt(lastResult.Qxx(3,3)) * C << " m" << endl;
        cout << "验后单位权中误差: " << sqrt(lastResult.sigma2) << endl;
    }
}
