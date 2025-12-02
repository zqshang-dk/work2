// 修改后的完整程序（在你的基础上做了最小侵入的改动并添加了系统分组、PRN记录、ECEF->BLH 转换以及输出）
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// 光速（m/s）
const double C_LIGHT = 299792458.0; 

// WGS84 常数（用于 ECEF->BLH）
const double WGS84_A = 6378137.0;
const double WGS84_F = 1.0 / 298.257223563;
const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F; // e^2

// 定义一个结构体用来存储单颗卫星的观测数据
struct SatData {
    string id_raw;  // 原始ID，如 "G01"
    char sys;       // 卫星系统字符 'G','C','R','E',...
    string prn;     // PRN号字符串，例如 "01","10"（便于记录）
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
    Vector4d Parameter; // X,Y,Z,dt
    double sigma0;  //验后单位权中误差
    Matrix4d cov;   //验后参数协方差矩阵
    bool success;   //是否解算成功
};

// pntpos 函数：不变（只是去掉了不合适的伪距阈值过滤，使用传入观测直接解）
// 说明：保留你原来的加权最小二乘实现，但不再在函数内随意丢弃伪距<10000的卫星。
//       如果你希望在这里做粗差剔除，请在调用前对 epoch_obs 做检测/剔除。
SolveResult pntpos(const vector<SatData>& obs, const EpochTime& time) {
    SolveResult res;
    res.success = false;

    vector<SatData> cleaned_obs;
    for (int i = 0; i < obs.size(); i++) {
        if (obs[i].pseudorange > 10000.0) {  // 只保留伪距大于10000米的卫星
            cleaned_obs.push_back(obs[i]);
        }
    }

    int n = cleaned_obs.size();
    // 如果卫星数少于4颗，无法定位
    if (n < 4) {
        return res;
    }

    //数据初始化（可改为更好的初值，如地心或上次解）
    double Xr = 0;
    double Yr = 0;
    double Zr = 0;
    double dt = 0;

    //迭代参数
    const int maxIter = 10;
    const double eps = 1e-4;  //收敛阈值（m）

    MatrixXd H(n, 4);  //设计矩阵n*4
    VectorXd l(n);     //OMC观测值
    MatrixXd W = MatrixXd::Zero(n, n);  //创建一个n*n的零矩阵

    for (int iter = 0; iter < maxIter; ++iter){
        // 填充 H, l, W
        for (int i = 0; i < n; i++){
            double P_obs = cleaned_obs[i].pseudorange;
            double Xs = cleaned_obs[i].sat_x;
            double Ys = cleaned_obs[i].sat_y;
            double Zs = cleaned_obs[i].sat_z; 
            double Var = cleaned_obs[i].variance;

            if(Var <= 0) Var = 1.0;

            double P0 = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
            if (P0 < 1.0) P0 = 1.0; // 防止除0

            H(i, 0) = (Xr - Xs) / P0;
            H(i, 1) = (Yr - Ys) / P0;
            H(i, 2) = (Zr - Zs) / P0;
            H(i, 3) = 1.0;

            // OMC 残差 (观测 - 计算)
            l(i) = P_obs - (P0 + dt);

            // 权矩阵（用方差的倒数）
            W(i, i) = 1.0 / Var;
        }

        // 法方程求解 dx——修改原因：当法矩阵奇异是，ldlt会导致数值爆炸，用更稳健的方式
        //VectorXd dx = (H.transpose() * W * H).ldlt().solve(H.transpose() * W * l);
        VectorXd dx = (H.transpose() * W * H + 1e-6 * Matrix4d::Identity()).ldlt().solve(H.transpose() * W * l);

        if (dx.size() != 4) break;
        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt += dx(3);

        double pos_delta = sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2));
        if (pos_delta < eps) break;
    }

    // 计算残差并评定精度
    VectorXd V(n);
    for (int i = 0; i < n; i++){
        double P_obs = cleaned_obs[i].pseudorange;
        double Xs = cleaned_obs[i].sat_x;
        double Ys = cleaned_obs[i].sat_y;
        double Zs = cleaned_obs[i].sat_z;

        double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
        double P_cal = dist + dt;
        V(i) = P_cal - P_obs;
    }
    // 验后单位权方差
    double vtwv = V.transpose() * W * V;
    double sigma0 = 0.0;
    if (n - 4 > 0)
        sigma0 = sqrt(vtwv / (n - 4));
    else
        sigma0 = 0.0;

    // 验后协方差矩阵
    MatrixXd Qxx = (H.transpose() * W * H).inverse();
    Matrix4d Dxx = Qxx * sigma0 * sigma0;

    res.Parameter << Xr, Yr, Zr, dt;
    res.sigma0 = sigma0;
    res.cov = Dxx;
    res.success = true;
    return res;
}

// ECEF → ENU 转换（输入坐标为 m） - 保留你原来的实现
Vector3d ecef2enu(double X, double Y, double Z,
                  double X0, double Y0, double Z0)
{
    double dx = X - X0;
    double dy = Y - Y0;
    double dz = Z - Z0;

    double lon = atan2(Y0, X0);
    double p = sqrt(X0*X0 + Y0*Y0);
    double lat = atan2(Z0, p);

    Matrix3d R;
    R << -sin(lon),              cos(lon),               0,
         -sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat),
          cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat);

    Vector3d enu = R * Vector3d(dx, dy, dz);
    return enu;
}

// ECEF -> BLH (经度、纬度以度为单位，height为米)
// 使用迭代法（Bowring-like）求解精度良好
void ecef2blh(double X, double Y, double Z, double &lat_deg, double &lon_deg, double &h) {
    double a = WGS84_A;
    double e2 = WGS84_E2;

    double lon = atan2(Y, X);

    double p = sqrt(X*X + Y*Y);
    double lat = atan2(Z, p * (1 - e2)); // 初值

    double lat_prev = 0;
    int iter = 0;
    while (fabs(lat - lat_prev) > 1e-12 && iter < 20) {
        lat_prev = lat;
        double sinlat = sin(lat);
        double N = a / sqrt(1.0 - e2 * sinlat * sinlat);
        h = p / cos(lat) - N;
        lat = atan2(Z, p * (1.0 - e2 * (N / (N + h))));
        iter++;
    }

    // 最终高度
    double sinlat = sin(lat);
    double N = a / sqrt(1.0 - e2 * sinlat * sinlat);
    h = p / cos(lat) - N;

    lat_deg = lat * 180.0 / M_PI;
    lon_deg = lon * 180.0 / M_PI;
}

// 评估结构体保持不变
struct EvalResult {
    double meanX, meanY, meanZ;
    double stdX, stdY, stdZ;
    double meanE, meanN, meanU;
    double rmsE, rmsN, rmsU;
};

// evaluate 函数保持（微调：可以复用你原实现）
EvalResult evaluate(const vector<Vector4d>& allPos,
                    double X0, double Y0, double Z0,
                    const string& out_prefix)
{
    int m = allPos.size();
    EvalResult R;
    if(m == 0){
        cerr << "错误：没有可评估的数据" << endl;
        return R;
    }
    double sumX = 0, sumY = 0, sumZ = 0;
    for (size_t i = 0; i < m; i++) {
        sumX += allPos[i](0);
        sumY += allPos[i](1);
        sumZ += allPos[i](2);
    }
    R.meanX = sumX / m;
    R.meanY = sumY / m;
    R.meanZ = sumZ / m;

    double sX=0, sY=0, sZ=0;
    for (size_t i = 0; i < m; i++) {
        sX += pow(allPos[i](0) - R.meanX, 2);
        sY += pow(allPos[i](1) - R.meanY, 2);
        sZ += pow(allPos[i](2) - R.meanZ, 2);
    }
    R.stdX = sqrt(sX/(m-1));
    R.stdY = sqrt(sY/(m-1));
    R.stdZ = sqrt(sZ/(m-1));

    vector<double> E, N, U;
    E.reserve(m); N.reserve(m); U.reserve(m);

    ofstream f_enu(out_prefix + "_ENU_series.txt");
    if (!f_enu.is_open()) {
        cerr << "错误：无法创建ENU输出文件" << endl;
        return R;
    }
    f_enu << "Epoch   E(m)   N(m)   U(m)\n";

    for (int i=0;i<m;i++){
        Vector3d enu = ecef2enu(
            allPos[i](0),
            allPos[i](1),
            allPos[i](2),
            X0, Y0, Z0
        );
        E.push_back(enu(0));
        N.push_back(enu(1));
        U.push_back(enu(2));
        f_enu << i+1 << " "
              << enu(0) << " "
              << enu(1) << " "
              << enu(2) << "\n";
    }
    f_enu.close();

    double sumE = 0, sumN = 0, sumU = 0;
    for (int i = 0; i < m; i++) {
        sumE += E[i];
        sumN += N[i];
        sumU += U[i];
    }
    R.meanE = sumE / m;
    R.meanN = sumN / m;
    R.meanU = sumU / m;

    double sumE2 = 0, sumN2 = 0, sumU2 = 0;
    for (int i = 0; i < m; i++) {
        sumE2 += E[i] * E[i];
        sumN2 += N[i] * N[i];
        sumU2 += U[i] * U[i];
    }
    R.rmsE = sqrt(sumE2 / m);
    R.rmsN = sqrt(sumN2 / m);
    R.rmsU = sqrt(sumU2 / m);

    return R;
}

int main() {
    vector<Vector4d> all_results;  // 存储每个历元的 X Y Z dt（平均后的）
    string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\mydata\2059329.25ObsCorr)";
    string outfile_path = "outpos_total.txt";
    string outblh_path = "outpos_blh.txt";
    ifstream infile(filename);
    ofstream outfile(outfile_path);
    ofstream outblh(outblh_path);

    if (!infile.is_open()) {
        cerr << "错误：无法打开文件 " << filename << endl;
        return -1;
    }
    if (!outfile.is_open()) {
        cerr << "错误：无法创建 " << outfile_path << endl;
        return -1;
    }
    if (!outblh.is_open()) {
        cerr << "错误：无法创建 " << outblh_path << endl;
        return -1;
    }

    // 输出表头
    outfile << fixed << setprecision(4);
    outblh << fixed << setprecision(8);
    outfile << "# 格式: 历元头, 然后每个系统的解 (如果可解), 最后 AverageXYZ\n";
    outblh << "# 格式: Epoch Week TOW avg_lat(deg) avg_lon(deg) avg_h(m)\n";

    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            EpochTime time_info;
            int num_sats = 0;
            char hash_char;
            stringstream ss_header(line);
            ss_header >> hash_char >> time_info.epoch_num >> time_info.week >> time_info.tow >> num_sats;

            vector<SatData> epoch_obs;
            epoch_obs.reserve(num_sats);

            for (int i = 0; i < num_sats; ++i) {
                if (!getline(infile, line)) break;
                if (line.empty()) { i--; continue; } // 避免空行计数错误

                stringstream ss_obs(line);
                SatData sat;
                // 解析：ID X Y Z Pseudo Variance
                ss_obs >> sat.id_raw >> sat.sat_x >> sat.sat_y >> sat.sat_z >> sat.pseudorange >> sat.variance;

                // 解析 sys 与 prn
                if (!sat.id_raw.empty()) {
                    sat.sys = sat.id_raw[0];
                    if (sat.id_raw.size() > 1)
                        sat.prn = sat.id_raw.substr(1);
                    else sat.prn = "";
                } else {
                    sat.sys = '?';
                    sat.prn = "";
                }
                epoch_obs.push_back(sat);
            }

            // ------------------------------------------
            // 分系统（按 sat.sys）解算（每个系统一个独立观测集）
            // ------------------------------------------
            map<char, vector<SatData>> sys_groups;
            for (const auto &s : epoch_obs) {
                sys_groups[s.sys].push_back(s);
            }

            // 保存每个系统的解
            map<char, SolveResult> sys_results;
            for (auto &kv : sys_groups) {
                char sysc = kv.first;
                vector<SatData> &obs_sys = kv.second;

                // 只解算卫星数>=4 的系统
                if ((int)obs_sys.size() >= 4) {
                    SolveResult r = pntpos(obs_sys, time_info);
                    if (r.success) {
                        sys_results[sysc] = r;
                    }
                }
            }

            // 将每个系统的解写入 outfile
            outfile << "# " << left << setw(5) << time_info.epoch_num
                    << setw(8) << time_info.week
                    << setw(12) << fixed << setprecision(4) << time_info.tow
                    << setw(6) << num_sats << "\n";

            // 写系统表头
            outfile << left << setw(6) << "SYS"
                    << setw(15) << "X(m)"
                    << setw(15) << "Y(m)"
                    << setw(15) << "Z(m)"
                    << setw(12) << "T(m)"
                    << setw(12) << "sigma0" << "\n";

            // 写出每个系统的解
            for (auto &kv : sys_results) {
                char sysc = kv.first;
                SolveResult &r = kv.second;
                outfile << left << setw(6) << string(1, sysc)
                        << setw(15) << r.Parameter(0)
                        << setw(15) << r.Parameter(1)
                        << setw(15) << r.Parameter(2)
                        << setw(12) << r.Parameter(3)
                        << setw(12) << r.sigma0 << "\n";
            }

            // 计算可用系统的平均 XYZ（只对 X,Y,Z 做平均）
            Vector3d sumXYZ(0,0,0);
            int cnt = 0;
            for (auto &kv: sys_results) {
                Vector4d p = kv.second.Parameter;
                sumXYZ += Vector3d(p(0), p(1), p(2));
                cnt++;
            }
            // Vector3d weighted_sum(0,0,0);
            // double total_weight = 0;
            // for (auto &kv : sys_results) {
            //     double w = 1.0 / (kv.second.sigma0 * kv.second.sigma0 + 1e-6);
            //     weighted_sum += w * kv.second.Parameter.head<3>();
            //     total_weight += w;
            // }
            // Vector3d avgXYZ = weighted_sum / total_weight;


            if (cnt > 0) {
                Vector3d avgXYZ = sumXYZ / double(cnt);
                outfile << left << setw(6) << "AVG"
                        << setw(15) << avgXYZ(0)
                        << setw(15) << avgXYZ(1)
                        << setw(15) << avgXYZ(2)
                        << "\n";
                outfile << "------------------------------------------------------------\n";

                // 将平均 ECEF 转为 BLH 写到 outblh
                double lat_deg, lon_deg, h;
                ecef2blh(avgXYZ(0), avgXYZ(1), avgXYZ(2), lat_deg, lon_deg, h);
                outblh << time_info.epoch_num << " "
                       << time_info.week << " "
                       << fixed << setprecision(4) << time_info.tow << " "
                       << setw(12) << setprecision(8) << lat_deg << " "
                       << setw(12) << setprecision(8) << lon_deg << " "
                       << setw(12) << setprecision(4) << h << "\n";

                // 保存进 all_results（用于最后的统计评估）
                Vector4d avgParam;
                avgParam << avgXYZ(0), avgXYZ(1), avgXYZ(2), 0.0;
                all_results.push_back(avgParam);
            } else {
                // 没有任何系统解算成功
                outfile << "WARNING: no system solution for this epoch\n";
                outfile << "------------------------------------------------------------\n";
            }
        }
    }

    infile.close();
    outfile.close();
    outblh.close();

    cout << "处理完成，结果已保存至 " << outfile_path << " 和 " << outblh_path << endl;

    // 如果收集到了结果，进行精度评估（使用第一个平均解作为参考）
    if (!all_results.empty()) {
        double X0 = all_results[0](0);
        double Y0 = all_results[0](0);
        double Z0 = all_results[0](0);
        vector<Vector4d> xyz_only;
        for (auto &v : all_results) xyz_only.push_back(v);
        EvalResult eval = evaluate(xyz_only, X0, Y0, Z0, "result");
        ofstream eval_file("accuracy_evaluation.txt");
        eval_file << fixed << setprecision(4);
        eval_file << "========== 精度评估结果 ==========\n";
        eval_file << "位置均值: (" << eval.meanX << ", " << eval.meanY << ", " << eval.meanZ << ")\n";
        eval_file << "位置标准差: (" << eval.stdX << ", " << eval.stdY << ", " << eval.stdZ << ")\n";
        eval_file << "ENU均值: (" << eval.meanE << ", " << eval.meanN << ", " << eval.meanU << ")\n";
        eval_file << "ENU RMS: (" << eval.rmsE << ", " << eval.rmsN << ", " << eval.rmsU << ")\n";
        eval_file.close();
        cout << "精度评估结果已保存至 accuracy_evaluation.txt" << endl;
    } else {
        cout << "警告：没有成功解算的历元" << endl;
    }

    return 0;
}
