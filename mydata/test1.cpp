#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//光速（m/s）
const double C_LIGHT = 299792458.0;
const double A = 6378137.0;
const double F = 1.0 / 298.257223563;
const double E2 = 2 * F - F * F;

// 定义一个结构体用来存储单颗卫星的观测数据
struct SatData {
    string id;      // 卫星ID (例如 C01)
    double sat_x;        // 卫星 X 坐标 (m)
    double sat_y;        // 卫星 Y 坐标 (m)
    double sat_z;        // 卫星 Z 坐标 (m)
    double pseudorange;  // 伪距观测值 (m) (对应文件第5列)
    double variance;     // 误差方差 (对应文件第6列)
};

// 定义一个结构体存储当前历元的时间信息
struct EpochTime {
    int epoch_num; // 历元编号
    int week;      // GPS周
    double tow;    // 周内秒 (Time of Week)
};

//定义解算结果结构体
struct SolveResult{
    VectorXd Parameter;
    double sigma0;  //验后单位权中误差
    //Matrix4d cov;   //验后参数协方差矩阵
    bool success;   //是否解算成功
};

// ECEF → BLH 转换函数（高精度迭代法）
void ecef2blh(double X, double Y, double Z, double& B, double& L, double& H) {
    double p = sqrt(X*X + Y*Y);
    L = atan2(Y, X);                    // 经度（rad）
    double B0 = atan2(Z, p * (1 - E2)); // 初始纬度
    double B_prev;
    double N;
    do {
        B_prev = B0;
        N = A / sqrt(1 - E2 * sin(B0)*sin(B0));
        H = p / cos(B0) - N;
        B0 = atan2(Z, p * (1 - E2 * N / (N + H)));
    } while (fabs(B0 - B_prev) > 1e-12);

    B = B0;
    N = A / sqrt(1 - E2 * sin(B)*sin(B));
    H = p / cos(B) - N;

    B = B * 180.0 / M_PI;   // 转为度
    L = L * 180.0 / M_PI;
}

SolveResult pntpos_multi(const vector<SatData>& obs) {
    SolveResult res;
    res.success = false;
    res.Parameter = VectorXd::Zero(7);          //// X,Y,Z, dtG, dtC, dtE, dtR
    vector<SatData> cleaned_obs;
    for (int i = 0; i < obs.size(); i++) {
        if (obs[i].pseudorange > 10000.0 && obs[i].variance>=0) {  // 只保留伪距大于10000米的卫星
            cleaned_obs.push_back(obs[i]);
        }
    }
    //太悲伤了！！！！！因为没有进行数据清洗导致02和03的粗差没有消除！
    
    int n = cleaned_obs.size();
    
    // 如果卫星数少于4颗，无法定位
    if (n < 7) {
        return res;
    }

    // VectorXd x = VectorXd::Zero(7);
    // x.head<3>() = Vector3d(0, 0, 0);

    
    //数据初始化
    double Xr = 0;
    double Yr = 0;
    double Zr = 0;
    double dt_GPS = 0;
    double dt_BDS = 0;
    double dt_GLO = 0;
    double dt_Gal = 0;

    //迭代参数
    const int maxIter = 15;
    const double eps = 1e-4;  //收敛阈值（m）
    
    MatrixXd H(n, 7);  //设计矩阵
    VectorXd l(n);     //OMC观测值
    MatrixXd W = MatrixXd::Zero(n, n);  //创建一个n*n的零矩阵

    for (int iter = 0; iter < maxIter; ++iter){
        //用Eigen来定义矩阵和向量    
        H.setZero();
        l.setZero();
        W.setZero();
        //填充H,l,W，此时需要遍历所有星历才可以实现
        for (int i = 0; i < n;i++){
            double P_obs = cleaned_obs[i].pseudorange;
            double Xs = cleaned_obs[i].sat_x;
            double Ys = cleaned_obs[i].sat_y;
            double Zs = cleaned_obs[i].sat_z; 
            double Var = cleaned_obs[i].variance;
            //double P0 = sqrt((Xs - x(0)) * (Xs - x(0)) + (Ys - x(1)) * (Ys - x(1)) + (Zs - x(2)) * (Zs - x(2)));
            
            
            if(Var<=0)
                Var = 1.0;
            
            double P0 = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
            //H.row(i).setZero();
            
            // H(i, 0) = (x(0) - Xs) / P0;
            // H(i, 1) = (x(1) - Ys) / P0;
            // H(i, 2) = (x(2) - Zs) / P0;
            H(i, 0) = (Xr - Xs) / P0;
            H(i, 1) = (Yr - Ys) / P0;
            H(i, 2) = (Zr - Zs) / P0;
            // H(i, 3) = 1.0;

            //OMC残差
            // l(i) = P_obs - (P0 + dt);

            //钟差列：根据系统设置1
            char sys = cleaned_obs[i].id[0];
            if(sys=='G')
                H(i, 3) = 1.0;
            else if(sys=='C')
                H(i, 4) = 1.0;
            else if(sys=='E')
                H(i, 5) = 1.0;
            else if(sys=='R')
                H(i, 6) = 1.0;
            else
                H(i, 3) = 1.0;  //如果没写的话默认是GPS卫星
            
            //观测方程中钟差的设置
            // double dt_sys = 0.0;
            // if (sys == 'G')
            //     dt_sys = x(3);
            // else if (sys == 'C')
            //     dt_sys = x(4);
            // else if (sys == 'E')
            //     dt_sys = x(5);
            // else if (sys == 'R')
            //     dt_sys = x(6);
            double dt_sys = 0.0;
            if (sys == 'G')
                dt_sys = dt_GPS;
            else if (sys == 'C')
                dt_sys = dt_BDS;
            else if (sys == 'E')
                dt_sys = dt_Gal;
            else if (sys == 'R')
                dt_sys = dt_GLO;

            //OMC残差
            l(i) = P_obs - (P0 + dt_sys);

            //构造权矩阵
            W(i, i) = 1.0 / Var;
        }
        VectorXd dx = (H.transpose() * W * H).ldlt() .solve (H.transpose() * W * l);

        // Xr += dx(0);
        // Yr += dx(1);
        // Zr += dx(2);
        // dt += dx(3);
        //x += dx;
        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt_GPS+=dx(3);
        dt_BDS+=dx(4);
        dt_Gal+=dx(5);
        dt_GLO += dx(6);

        double pos_delta = sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2));
        if(pos_delta<eps){
            break;
        }
    }

    //精度评定
    VectorXd V(n);  //残差
    for (int i = 0; i < n;i++){
        double P_obs = cleaned_obs[i].pseudorange;
        double Xs = cleaned_obs[i].sat_x;
        double Ys = cleaned_obs[i].sat_y;
        double Zs = cleaned_obs[i].sat_z;
        char sys = cleaned_obs[i].id[0];
        // double dt_sys = (sys=='G'?x(3):(sys=='C'?x(4):(sys=='E'?x(5):x(6))));
        // double dist = sqrt((Xs - x(0)) * (Xs - x(0)) + (Ys - x(1)) * (Ys - x(1)) + (Zs - x(2)) * (Zs - x(2)));
        double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
        double dt_sys = (sys=='G'?dt_GPS:(sys=='C'?dt_BDS:(sys=='E'?dt_Gal:dt_GLO)));
        double P_cal = dist + dt_sys;
        V(i) = P_cal - P_obs;
    }
    //计算验后单位权方差
    double vtwv = V.transpose() * W * V;
    double sigma0 = sqrt(vtwv / (n - 7));

    // //计算验后协方差矩阵
    // MatrixXd Qxx= (H.transpose() * W * H).inverse();

    // Matrix4d Dxx = Qxx * sigma0 * sigma0;

    //res.Parameter << x(0), x(1), x(2), dt;
    // res.Parameter.head<3>() = x.head<3>();
    // res.Parameter.tail<4>() = x.tail<4>();
    res.Parameter << Xr, Yr, Zr, dt_GPS, dt_BDS, dt_Gal, dt_GLO;
    res.sigma0 = sigma0;
    //res.cov = Dxx;
    res.success = true;

    return res;
}

// // ECEF → ENU 转换（输入坐标为 m）
// Vector3d ecef2enu(double X, double Y, double Z,
//                   double X0, double Y0, double Z0)
// {
//     // 计算参考点的经纬度
//     double dx = X - X0;
//     double dy = Y - Y0;
//     double dz = Z - Z0;

//     double lon = atan2(Y0, X0);
//     double p = sqrt(X0*X0 + Y0*Y0);
//     double lat = atan2(Z0, p);

//     // 构建旋转矩阵
//     Matrix3d R;
//     R << -sin(lon),              cos(lon),               0,
//          -sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat),
//           cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat);

//     Vector3d enu = R * Vector3d(dx, dy, dz);
//     return enu;
// }

struct EvalResult {
    double meanX, meanY, meanZ;
    double stdX, stdY, stdZ;
    double meanE, meanN, meanU;
    double rmsE, rmsN, rmsU;
};

EvalResult evaluate(const vector<VectorXd>& allPos,
                    double X0, double Y0, double Z0,
                    const string& out_prefix)
{
    int m = allPos.size();
    EvalResult R;

    if(m==0){
        cerr << "错误：没有可评估的数据" << endl;
        return R;
    }
    // ================================
    // 1) 计算 XYZ 均值
    // ================================
    double sumX = 0, sumY = 0, sumZ = 0;
    for (size_t i = 0; i < m; i++) {
        sumX += allPos[i](0);
        sumY += allPos[i](1);
        sumZ += allPos[i](2);
    }
    R.meanX = sumX / m;
    R.meanY = sumY / m;
    R.meanZ = sumZ / m;

    // ================================
    // 2) 样本标准差
    // ================================
    double sX=0, sY=0, sZ=0;
    for (size_t i = 0; i < m; i++) {
        sX += pow(allPos[i](0) - R.meanX, 2);
        sY += pow(allPos[i](1) - R.meanY, 2);
        sZ += pow(allPos[i](2) - R.meanZ, 2);
    }
    R.stdX = sqrt(sX/(m-1));
    R.stdY = sqrt(sY/(m-1));
    R.stdZ = sqrt(sZ/(m-1));

    // // ================================
    // // 3) 计算每个历元的 ENU
    // // ================================
    // vector<double> E, N, U;
    // E.reserve(m); N.reserve(m); U.reserve(m);

    // ofstream f_enu(out_prefix + "_ENU_series.txt");
    // if (!f_enu.is_open()) {
    //     cerr << "错误：无法创建ENU输出文件" << endl;
    //     return R;
    // }

    // f_enu << "Epoch   E(m)   N(m)   U(m)\n";

    // for (int i=0;i<m;i++){
    //     Vector3d enu = ecef2enu(
    //         allPos[i](0),
    //         allPos[i](1),
    //         allPos[i](2),
    //         X0, Y0, Z0
    //     );

    //     E.push_back(enu(0));
    //     N.push_back(enu(1));
    //     U.push_back(enu(2));

    //     f_enu << i+1 << " "
    //           << enu(0) << " "
    //           << enu(1) << " "
    //           << enu(2) << "\n";
    // }
    // f_enu.close();

    // ================================
    // 4) ENU 均值 + RMS
    // ================================
    // auto mean = [&](vector<double>& v){
    //     double s=0; for(double x:v) s+=x; return s/m;
    // };
    // auto rms = [&](vector<double>& v){
    //     double s=0; for(double x:v) s+=x*x; return sqrt(s/m);
    // };

    // R.meanE = mean(E);
    // R.meanN = mean(N);
    // R.meanU = mean(U);

    // R.rmsE = rms(E);
    // R.rmsN = rms(N);
    // R.rmsU = rms(U);
    // double sumE = 0, sumN = 0, sumU = 0;
    // for (int i = 0; i < m; i++) {
    //     sumE += E[i];
    //     sumN += N[i];
    //     sumU += U[i];
    // }
    // R.meanE = sumE / m;
    // R.meanN = sumN / m;
    // R.meanU = sumU / m;
    
    // // 计算RMS
    // double sumE2 = 0, sumN2 = 0, sumU2 = 0;
    // for (int i = 0; i < m; i++) {
    //     sumE2 += E[i] * E[i];
    //     sumN2 += N[i] * N[i];
    //     sumU2 += U[i] * U[i];
    // }
    // R.rmsE = sqrt(sumE2 / m);
    // R.rmsN = sqrt(sumN2 / m);
    // R.rmsU = sqrt(sumU2 / m);

    // return R;
}


int main() {
    vector<VectorXd> all_results;  // 存储每个历元的 X Y Z dt

    // 原始文件路径 (使用了 R"()" 语法，不需要双斜杠)
    string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\mydata\2059330.25ObsCorr)";
    string outfile_path = "outpos_total.txt";
    string blhfile_path = "blh_results.txt";

    ifstream infile(filename);
    ofstream outfile(outfile_path);
    ofstream blhfile(blhfile_path);  // 新增

    if (!infile.is_open()) {
        cerr << "错误：无法打开文件 " << filename << endl;
        return -1;
    }
     
    if (!blhfile.is_open()) {
        cerr << "错误：无法创建BLH输出文件" << endl;
        return -1;
    }

    // 设置输出精度为小数点后3位
    outfile << fixed << setprecision(3);
    blhfile << fixed << setprecision(8); 

    blhfile << "# 格式: Epoch Week TOW avg_lat(deg) avg_lon(deg) avg_h(m)\n";

    string line;
    while (getline(infile, line)) {
        // 跳过空行
        if (line.empty()) continue;

        // 检查是否是历元头 (以 '#' 开头)
        if (line[0] == '#') {
            EpochTime time_info;
            int num_sats = 0;
            char hash_char; 

            // 解析头文件: "#  1 2170 172800.000 25"
            stringstream ss_header(line);
            ss_header >> hash_char >> time_info.epoch_num >> time_info.week >> time_info.tow >> num_sats;

            // 读取该历元下的所有卫星数据
            vector<SatData> epoch_obs;
            epoch_obs.reserve(num_sats); 

            for (int i = 0; i < num_sats; ++i) {
                if (!getline(infile, line)) break;
                
                stringstream ss_obs(line);
                SatData sat;
                // 解析每行卫星数据
                ss_obs >> sat.id >> sat.sat_x >> sat.sat_y >> sat.sat_z >> sat.pseudorange >> sat.variance;
                
                epoch_obs.push_back(sat);
            }

            // 进行定位解算
            if (epoch_obs.size() >= 4) {
                vector<double> result(4, 0.0); // 存放解 [x, y, z, dt]
                
                SolveResult res = pntpos_multi(epoch_obs);

                vector<VectorXd> all_positions;

                if (res.success) {
                    //all_results.push_back(res.Parameter);

                    double X = res.Parameter(0);
                    double Y = res.Parameter(1);
                    double Z = res.Parameter(2);
                    double Lat, Lon, H;
                    ecef2blh(X, Y, Z, Lat, Lon, H);

                    // outfile << "# " << time_info.epoch_num << " "
                    //         << time_info.week << " " << fixed << setprecision(3) << time_info.tow << " "
                    //         << epoch_obs.size() << " sats" << endl;

                    // outfile << fixed << setprecision(4);
                    // outfile << "X: " << setw(15) << X
                    //         << " Y: " << setw(15) << Y
                    //         << " Z: " << setw(15) << Z << endl;

                    // outfile << "Lat: " << setw(12) << setprecision(8) << Lat << " deg"
                    //         << " Lon: " << setw(12) << Lon << " deg"
                    //         << " H: " << setw(10) << setprecision(3) << H << " m"
                    //         << "  sigma0: " << res.sigma0 << " m" << endl;

                    // outfile << "dt_GPS: " << setw(10) << setprecision(3) << res.Parameter(3)
                    //         << " dt_BDS: " << res.Parameter(4)
                    //         << " dt_GAL: " << res.Parameter(5)
                    //         << " dt_GLO: " << res.Parameter(6) << " m" << endl;
                    // ===== 输出头部 =====
                    outfile << "# " << setw(4) << time_info.epoch_num
                            << setw(8) << time_info.week
                            << setw(14) << fixed << setprecision(4) << time_info.tow
                            << setw(4) << epoch_obs.size() << "\n";

                    // ===== 输出标题行 =====
                    outfile << "X (m)          Y (m)          Z (m)          "
                            << "dT_GPS(m)      dT_BDS(m)      dT_GAL(m)      dT_GLO(m)\n";

                    // ===== 输出数据行 =====
                    outfile << fixed << setprecision(4)
                            << setw(15) << X
                            << setw(15) << Y
                            << setw(15) << Z
                            << setw(15) << res.Parameter(3)   // dT_GPS
                            << setw(15) << res.Parameter(4)   // dT_BDS
                            << setw(15) << res.Parameter(5)   // dT_GAL
                            << setw(15) << res.Parameter(6)   // dT_GLO
                            << "\n\n";

                    outfile << "------------------------------------------------------------" << endl;
                    // // 第一行：历元编号 GPS周 周秒 观测值数
                    // outfile << "# " << left << setw(5) << time_info.epoch_num 
                    //         << setw(8) << time_info.week 
                    //         << setw(12) << fixed << setprecision(4) << time_info.tow 
                    //         << num_sats << endl;

                    // // 第二行：表头
                    // outfile << left << setw(15) << "X (m)" 
                    //         << setw(15) << "Y (m)" 
                    //         << setw(15) << "Z (m)" 
                    //         << "T(m)" << endl;

                    // // 第三行：解算数值
                    // outfile << left << fixed << setprecision(4) 
                    //         << setw(15) << res.Parameter(0) 
                    //         << setw(15) << res.Parameter(1) 
                    //         << setw(15) << res.Parameter(2) 
                    //         << res.Parameter(3) << endl;

                    // // 第四行：验后单位权中误差
                    // outfile << "验后单位权中误差: " << res.sigma0 << " (m)" << endl;

                    // // 第五行：验后估计方差标题
                    // outfile << "验后估计方差(m^2)" << endl;

                    // // 第六至九行：4x4 协方差矩阵
                    // for (int r = 0; r < 4; ++r) {
                    //     outfile << left << fixed << setprecision(4)
                    //             << setw(15) << res.cov(r, 0)
                    //             << setw(15) << res.cov(r, 1)
                    //             << setw(15) << res.cov(r, 2)
                    //             << setw(15) << res.cov(r, 3) << endl;
                    // }

                    // // 分隔线
                    // outfile << "------------------------------------------------------------" << endl;

                    // ===== 新增：输出到独立的BLH文件 =====
                    blhfile << setw(5) << time_info.epoch_num  // Epoch
                            << setw(6) << time_info.week       // Week
                            << setw(14) << fixed << setprecision(4) << time_info.tow  // TOW
                            << "  "  // 两个空格分隔
                            << setw(15) << fixed << setprecision(8) << Lat    // Latitude
                            << " "
                            << setw(15) << fixed << setprecision(8) << Lon    // Longitude
                            << " "
                            << setw(12) << fixed << setprecision(4) << H      // Height
                            << "\n";

                    // 保存历元结果
                    all_results.push_back(res.Parameter);
                }
            }
        }
    }


    infile.close();
    outfile.close();
    blhfile.close();
    cout << "处理完成，结果已保存至 " << outfile_path << endl;

    // 如果收集到了结果，进行精度评估
    if (!all_results.empty()) {
        // 设置参考坐标（可以设为第一个历元的解或真实坐标）
        double X0 = -1132914.6126;  // 从你的第一个历元结果
        double Y0 = 6092528.9442;
        double Z0 = 1504633.0098;
        
        EvalResult eval = evaluate(all_results, X0, Y0, Z0, "result");
        
        // 输出评估结果
        // ofstream eval_file("accuracy_evaluation.txt");
        // eval_file << fixed << setprecision(4);
        // eval_file << "========== 精度评估结果 ==========\n";
        // eval_file << "位置均值: (" << eval.meanX << ", " << eval.meanY << ", " << eval.meanZ << ")\n";
        // eval_file << "位置标准差: (" << eval.stdX << ", " << eval.stdY << ", " << eval.stdZ << ")\n";
        // // eval_file << "ENU均值: (" << eval.meanE << ", " << eval.meanN << ", " << eval.meanU << ")\n";
        // // eval_file << "ENU RMS: (" << eval.rmsE << ", " << eval.rmsN << ", " << eval.rmsU << ")\n";
        // eval_file.close();
        
        cout << "精度评估结果已保存至 accuracy_evaluation.txt" << endl;
    } else {
        cout << "警告：没有成功解算的历元" << endl;
    }

    return 0;
}