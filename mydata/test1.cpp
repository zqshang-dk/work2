// 修改后的完整程序（修复E系统数据解析问题）
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

// WGS84 常数
const double WGS84_A = 6378137.0;
const double WGS84_F = 1.0 / 298.257223563;
const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F;

// 卫星观测数据结构
struct SatData {
    string id_raw;
    char sys;
    string prn;
    double sat_x, sat_y, sat_z;
    double pseudorange;
    double variance;
};

// 历元时间结构
struct EpochTime {
    int epoch_num;
    int week;
    double tow;
};

// 解算结果结构
struct SolveResult{
    Vector4d Parameter;
    double sigma0;
    Matrix4d cov;
    bool success;
};

// 坐标合理性检查函数
bool isReasonableCoordinate(double x, double y, double z) {
    // 地球坐标合理范围检查（米）
    const double MAX_VALID = 100000000.0;  // 1亿米
    const double MIN_VALID = -100000000.0;
    
    return (x > MIN_VALID && x < MAX_VALID &&
            y > MIN_VALID && y < MAX_VALID &&
            z > MIN_VALID && z < MAX_VALID);
}

// 伪距合理性检查
bool isReasonablePseudorange(double pr) {
    return (pr > 1000000.0 && pr < 100000000.0);  // 1000km 到 10万km
}

// 改进的pntpos函数，增加数据验证
SolveResult pntpos(const vector<SatData>& obs, const EpochTime& time) {
    SolveResult res;
    res.success = false;

    int n = obs.size();
    if (n < 4) {
        return res;
    }

    // 数据验证：检查每个卫星的坐标和伪距是否合理
    vector<SatData> valid_obs;
    for (const auto& sat : obs) {
        if (isReasonableCoordinate(sat.sat_x, sat.sat_y, sat.sat_z) && 
            isReasonablePseudorange(sat.pseudorange)) {
            valid_obs.push_back(sat);
        }
    }
    
    if (valid_obs.size() < 4) {
        cerr << "历元 " << time.epoch_num << ": 有效卫星数不足 (" << valid_obs.size() << ")" << endl;
        return res;
    }

    double Xr = 0;
    double Yr = 0;
    double Zr = 0;
    double dt = 0;

    const int maxIter = 10;
    const double eps = 1e-4;

    int valid_n = valid_obs.size();
    MatrixXd H(valid_n, 4);
    VectorXd l(valid_n);
    MatrixXd W = MatrixXd::Zero(valid_n, valid_n);

    for (int iter = 0; iter < maxIter; ++iter){
        for (int i = 0; i < valid_n; i++){
            double P_obs = valid_obs[i].pseudorange;
            double Xs = valid_obs[i].sat_x;
            double Ys = valid_obs[i].sat_y;
            double Zs = valid_obs[i].sat_z; 
            double Var = valid_obs[i].variance;

            if(Var <= 0) Var = 1.0;

            double P0 = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
            if (P0 < 1.0) P0 = 1.0;

            H(i, 0) = (Xr - Xs) / P0;
            H(i, 1) = (Yr - Ys) / P0;
            H(i, 2) = (Zr - Zs) / P0;
            H(i, 3) = 1.0;

            l(i) = P_obs - (P0 + dt);
            W(i, i) = 1.0 / Var;
        }

        VectorXd dx = (H.transpose() * W * H).ldlt().solve(H.transpose() * W * l);

        if (dx.size() != 4) break;
        
        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt += dx(3);

        double pos_delta = sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2));
        if (pos_delta < eps) break;
    }

    // 计算残差
    VectorXd V(valid_n);
    for (int i = 0; i < valid_n; i++){
        double P_obs = valid_obs[i].pseudorange;
        double Xs = valid_obs[i].sat_x;
        double Ys = valid_obs[i].sat_y;
        double Zs = valid_obs[i].sat_z;

        double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
        double P_cal = dist + dt;
        V(i) = P_cal - P_obs;
    }

    double vtwv = V.transpose() * W * V;
    double sigma0 = (valid_n - 4 > 0) ? sqrt(vtwv / (valid_n - 4)) : 0.0;

    MatrixXd Qxx = (H.transpose() * W * H).inverse();
    Matrix4d Dxx = Qxx * sigma0 * sigma0;

    // 最终结果验证
    if (!isReasonableCoordinate(Xr, Yr, Zr)) {
        cerr << "历元 " << time.epoch_num << ": 解算坐标不合理" << endl;
        return res;
    }

    res.Parameter << Xr, Yr, Zr, dt;
    res.sigma0 = sigma0;
    res.cov = Dxx;
    res.success = true;
    return res;
}

// ECEF -> BLH 转换（保持不变）
void ecef2blh(double X, double Y, double Z, double &lat_deg, double &lon_deg, double &h) {
    double a = WGS84_A;
    double e2 = WGS84_E2;

    double lon = atan2(Y, X);
    double p = sqrt(X*X + Y*Y);
    double lat = atan2(Z, p * (1 - e2));

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

    double sinlat = sin(lat);
    double N = a / sqrt(1.0 - e2 * sinlat * sinlat);
    h = p / cos(lat) - N;

    lat_deg = lat * 180.0 / M_PI;
    lon_deg = lon * 180.0 / M_PI;
}

int main() {
    vector<Vector4d> all_results;
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
    if (!outfile.is_open() || !outblh.is_open()) {
        cerr << "错误：无法创建输出文件" << endl;
        return -1;
    }

    outfile << fixed << setprecision(4);
    outblh << fixed << setprecision(8);
    outfile << "# 格式: 历元头, 然后每个系统的解, 最后 AverageXYZ\n";
    outblh << "# 格式: Epoch Week TOW avg_lat(deg) avg_lon(deg) avg_h(m)\n";

    string line;
    int total_epochs = 0;
    int solved_epochs = 0;

    while (getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            EpochTime time_info;
            int num_sats = 0;
            char hash_char;
            stringstream ss_header(line);
            ss_header >> hash_char >> time_info.epoch_num >> time_info.week >> time_info.tow >> num_sats;

            total_epochs++;

            vector<SatData> epoch_obs;
            epoch_obs.reserve(num_sats);

            for (int i = 0; i < num_sats; ++i) {
                if (!getline(infile, line)) break;
                if (line.empty()) { i--; continue; }

                stringstream ss_obs(line);
                SatData sat;
                ss_obs >> sat.id_raw >> sat.sat_x >> sat.sat_y >> sat.sat_z >> sat.pseudorange >> sat.variance;

                // 数据验证和记录
                if (!isReasonableCoordinate(sat.sat_x, sat.sat_y, sat.sat_z)) {
                    cerr << "历元 " << time_info.epoch_num << " 卫星 " << sat.id_raw 
                         << " 坐标异常: (" << sat.sat_x << ", " << sat.sat_y << ", " << sat.sat_z << ")" << endl;
                    continue;  // 跳过这个异常卫星
                }

                if (!sat.id_raw.empty()) {
                    sat.sys = sat.id_raw[0];
                    sat.prn = (sat.id_raw.size() > 1) ? sat.id_raw.substr(1) : "";
                } else {
                    sat.sys = '?';
                    sat.prn = "";
                }
                epoch_obs.push_back(sat);
            }

            // 分系统解算
            map<char, vector<SatData>> sys_groups;
            for (const auto &s : epoch_obs) {
                sys_groups[s.sys].push_back(s);
            }

            map<char, SolveResult> sys_results;
            for (auto &kv : sys_groups) {
                char sysc = kv.first;
                vector<SatData> &obs_sys = kv.second;

                if (obs_sys.size() >= 4) {
                    SolveResult r = pntpos(obs_sys, time_info);
                    if (r.success) {
                        sys_results[sysc] = r;
                    }
                }
            }

            outfile << "# " << left << setw(5) << time_info.epoch_num
                    << setw(8) << time_info.week
                    << setw(12) << time_info.tow
                    << setw(6) << num_sats << "\n";

            outfile << left << setw(6) << "SYS"
                    << setw(15) << "X(m)" << setw(15) << "Y(m)" << setw(15) << "Z(m)"
                    << setw(12) << "T(m)" << setw(12) << "sigma0" << "\n";

            // 计算平均坐标（只使用有效的系统解）
            Vector3d sumXYZ(0,0,0);
            int valid_sys_count = 0;
            
            for (auto &kv : sys_results) {
                char sysc = kv.first;
                SolveResult &r = kv.second;
                
                outfile << left << setw(6) << string(1, sysc)
                        << setw(15) << r.Parameter(0) << setw(15) << r.Parameter(1) 
                        << setw(15) << r.Parameter(2) << setw(12) << r.Parameter(3)
                        << setw(12) << r.sigma0 << "\n";
                
                // 只累加合理的坐标
                if (isReasonableCoordinate(r.Parameter(0), r.Parameter(1), r.Parameter(2))) {
                    sumXYZ += Vector3d(r.Parameter(0), r.Parameter(1), r.Parameter(2));
                    valid_sys_count++;
                }
            }

            if (valid_sys_count > 0) {
                Vector3d avgXYZ = sumXYZ / double(valid_sys_count);
                outfile << left << setw(6) << "AVG"
                        << setw(15) << avgXYZ(0) << setw(15) << avgXYZ(1) 
                        << setw(15) << avgXYZ(2) << "\n";
                
                // 转换为BLH坐标
                double lat_deg, lon_deg, h;
                ecef2blh(avgXYZ(0), avgXYZ(1), avgXYZ(2), lat_deg, lon_deg, h);
                outblh << time_info.epoch_num << " " << time_info.week << " " 
                       << time_info.tow << " " << lat_deg << " " << lon_deg << " " << h << "\n";

                Vector4d avgParam;
                avgParam << avgXYZ(0), avgXYZ(1), avgXYZ(2), 0.0;
                all_results.push_back(avgParam);
                solved_epochs++;
            } else {
                outfile << "WARNING: 无有效系统解\n";
            }
            outfile << "------------------------------------------------------------\n";
        }
    }

    infile.close();
    outfile.close();
    outblh.close();

    cout << "处理完成！总历元: " << total_epochs << ", 成功解算: " << solved_epochs << endl;
    cout << "结果文件: " << outfile_path << " 和 " << outblh_path << endl;

    return 0;
}