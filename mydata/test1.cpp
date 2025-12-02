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

// 改进的pntpos函数 - 添加高度角权重和稳健估计
SolveResult improved_pntpos(const vector<SatData>& obs, const EpochTime& time) {
    SolveResult res;
    res.success = false;

    // 数据预处理：基本筛选
    vector<SatData> valid_obs;
    for (const auto& sat : obs) {
        if (sat.pseudorange < 1e6 || sat.pseudorange > 1e8) continue;
        if (!isReasonableCoordinate(sat.sat_x, sat.sat_y, sat.sat_z)) continue;
        valid_obs.push_back(sat);
    }
    
    int n = valid_obs.size();
    if (n < 4) {
        cerr << "历元 " << time.epoch_num << ": 有效卫星数不足 (" << n << ")" << endl;
        return res;
    }

    // 初始位置：使用所有卫星的几何中心
    double Xr = 0, Yr = 0, Zr = 0;
    for (const auto& sat : valid_obs) {
        Xr += sat.sat_x;
        Yr += sat.sat_y;
        Zr += sat.sat_z;
    }
    Xr /= n; Yr /= n; Zr /= n;
    double dt = 0;

    const int maxIter = 15;
    const double eps = 1e-4;
    
    Vector4d best_solution;
    double best_sigma0 = 1e9;
    int best_iteration = -1;

    // 迭代解算
    for (int iter = 0; iter < maxIter; ++iter) {
        MatrixXd H(n, 4);
        VectorXd l(n);
        VectorXd weights(n);
        
        // 计算设计矩阵和残差
        for (int i = 0; i < n; i++) {
            const auto& sat = valid_obs[i];
            double dx = Xr - sat.sat_x;
            double dy = Yr - sat.sat_y;
            double dz = Zr - sat.sat_z;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < 1.0) dist = 1.0;
            
            H(i, 0) = dx / dist;
            H(i, 1) = dy / dist;
            H(i, 2) = dz / dist;
            H(i, 3) = 1.0;
            
            l(i) = sat.pseudorange - (dist + dt);
            
            // 改进的权重计算：考虑高度角和信噪比
            double elevation = calculateElevation(Xr, Yr, Zr, sat.sat_x, sat.sat_y, sat.sat_z);
            double elevation_weight = sin(max(elevation, 5.0 * M_PI/180.0)); // 最小5度
            
            // 方差权重 + 高度角权重
            double variance = max(sat.variance, 1.0);
            weights(i) = elevation_weight * elevation_weight / variance;
        }
        
        // 粗差探测：基于IQR方法
        VectorXd residuals = l;
        if (iter > 0) { // 从第二次迭代开始粗差探测
            vector<double> sorted_residuals(residuals.data(), residuals.data() + n);
            sort(sorted_residuals.begin(), sorted_residuals.end());
            
            double Q1 = sorted_residuals[n/4];
            double Q3 = sorted_residuals[3*n/4];
            double IQR = Q3 - Q1;
            double lower_bound = Q1 - 2.0 * IQR;  // 宽松的阈值
            double upper_bound = Q3 + 2.0 * IQR;
            
            for (int i = 0; i < n; i++) {
                if (residuals(i) < lower_bound || residuals(i) > upper_bound) {
                    weights(i) *= 0.1; // 大幅降低粗差权重
                }
            }
        }
        
        // 加权最小二乘
        MatrixXd W = weights.asDiagonal();
        Matrix4d N = H.transpose() * W * H;
        
        // 检查法矩阵病态性
        Eigen::JacobiSVD<Matrix4d> svd(N);
        double cond = svd.singularValues()(0) / svd.singularValues()(3);
        if (cond > 1e6) {
            cerr << "历元 " << time.epoch_num << ": 法方程病态，条件数=" << cond << endl;
            break;
        }
        
        Vector4d dx = N.ldlt().solve(H.transpose() * W * l);
        
        if (!dx.allFinite()) {
            cerr << "历元 " << time.epoch_num << ": 解算出现无穷大值" << endl;
            break;
        }
        
        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt += dx(3);
        
        // 计算当前解的精度
        double sigma0 = calculateSigma0(valid_obs, Xr, Yr, Zr, dt, weights, n);
        
        // 记录最佳解
        if (sigma0 < best_sigma0 && sigma0 > 0.1) {
            best_sigma0 = sigma0;
            best_solution << Xr, Yr, Zr, dt;
            best_iteration = iter;
        }
        
        double pos_delta = dx.head<3>().norm();
        if (pos_delta < eps) {
            break;
        }
    }
    
    if (best_iteration >= 0) {
        res.Parameter = best_solution;
        res.sigma0 = best_sigma0;
        
        // 计算最终协方差矩阵
        MatrixXd H_final(n, 4);
        VectorXd weights_final(n);
        for (int i = 0; i < n; i++) {
            const auto& sat = valid_obs[i];
            double dx = best_solution(0) - sat.sat_x;
            double dy = best_solution(1) - sat.sat_y;
            double dz = best_solution(2) - sat.sat_z;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            H_final(i, 0) = dx / dist;
            H_final(i, 1) = dy / dist;
            H_final(i, 2) = dz / dist;
            H_final(i, 3) = 1.0;
            
            double elevation = calculateElevation(best_solution(0), best_solution(1), 
                                                best_solution(2), sat.sat_x, sat.sat_y, sat.sat_z);
            double elevation_weight = sin(max(elevation, 5.0 * M_PI/180.0));
            double variance = max(sat.variance, 1.0);
            weights_final(i) = elevation_weight * elevation_weight / variance;
        }
        
        MatrixXd W_final = weights_final.asDiagonal();
        Matrix4d Qxx = (H_final.transpose() * W_final * H_final).inverse();
        res.cov = Qxx * best_sigma0 * best_sigma0;
        res.success = true;
        
        cout << "历元 " << time.epoch_num << ": 最优解 sigma0=" << best_sigma0 
             << "m, 迭代=" << best_iteration << endl;
    }
    
    return res;
}

// 计算高度角辅助函数
double calculateElevation(double Xr, double Yr, double Zr, 
                         double Xs, double Ys, double Zs) {
    double dx = Xs - Xr;
    double dy = Ys - Yr;
    double dz = Zs - Zr;
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    // 计算接收机位置的大地坐标
    double lat_deg, lon_deg, h;
    ecef2blh(Xr, Yr, Zr, lat_deg, lon_deg, h);
    double lat = lat_deg * M_PI / 180.0;
    double lon = lon_deg * M_PI / 180.0;
    
    // 计算卫星在接收机本地坐标系中的向量
    double sin_lat = sin(lat), cos_lat = cos(lat);
    double sin_lon = sin(lon), cos_lon = cos(lon);
    
    double e = -sin_lon * dx + cos_lon * dy;
    double n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    double u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
    
    double elevation = atan2(u, sqrt(e*e + n*n));
    return elevation;
}

// 计算验后单位权中误差
double calculateSigma0(const vector<SatData>& obs, double Xr, double Yr, double Zr, double dt, 
                      const VectorXd& weights, int n) {
    VectorXd V(n);
    for (int i = 0; i < n; i++) {
        const auto& sat = obs[i];
        double dx = Xr - sat.sat_x;
        double dy = Yr - sat.sat_y;
        double dz = Zr - sat.sat_z;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        double residual = sat.pseudorange - (dist + dt);
        V(i) = residual;
    }
    
    double vtwv = V.transpose() * weights.asDiagonal() * V;
    return (n > 4) ? sqrt(vtwv / (n - 4)) : 0.0;
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
                    SolveResult r = improved_pntpos(obs_sys, time_info);

                    if (r.success) {
                        // 检查DOP值
                        MatrixXd H_local(obs_sys.size(), 4);
                        for (int i = 0; i < obs_sys.size(); i++) {
                            const auto& sat = obs_sys[i];
                            double dx = r.Parameter(0) - sat.sat_x;
                            double dy = r.Parameter(1) - sat.sat_y;
                            double dz = r.Parameter(2) - sat.sat_z;
                            double dist = sqrt(dx*dx + dy*dy + dz*dz);
                            
                            H_local(i, 0) = dx / dist;
                            H_local(i, 1) = dy / dist;
                            H_local(i, 2) = dz / dist;
                            H_local(i, 3) = 1.0;
                        }
                        
                        Matrix4d Q = (H_local.transpose() * H_local).inverse();
                        double pdop = sqrt(Q(0,0) + Q(1,1) + Q(2,2));
                        
                        if (pdop > 10.0 || r.sigma0 > 20.0) {
                            cout << "历元 " << time_info.epoch_num << " 系统 " << sysc 
                                << ": 解算质量较差 (PDOP=" << pdop << ", sigma0=" << r.sigma0 << "m)" << endl;
                            // 可以选择标记为不可用
                            // r.success = false;
                        }
                    }

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

            // // 计算平均坐标（只使用有效的系统解）
            // Vector3d sumXYZ(0,0,0);
            // int valid_sys_count = 0;
            
            // for (auto &kv : sys_results) {
            //     char sysc = kv.first;
            //     SolveResult &r = kv.second;
                
            //     outfile << left << setw(6) << string(1, sysc)
            //             << setw(15) << r.Parameter(0) << setw(15) << r.Parameter(1) 
            //             << setw(15) << r.Parameter(2) << setw(12) << r.Parameter(3)
            //             << setw(12) << r.sigma0 << "\n";
                
            //     // 只累加合理的坐标
            //     if (isReasonableCoordinate(r.Parameter(0), r.Parameter(1), r.Parameter(2))) {
            //         sumXYZ += Vector3d(r.Parameter(0), r.Parameter(1), r.Parameter(2));
            //         valid_sys_count++;
            //     }
            // }
            // ★ 第一步：先进行系统间一致性检查 ★
            vector<Vector3d> valid_solutions;
            vector<char> valid_systems;
            vector<double> valid_sigma0;

            for (auto &kv : sys_results) {
                SolveResult &r = kv.second;
                if (r.success && r.sigma0 < 15.0) { // 质量筛选
                    valid_solutions.push_back(r.Parameter.head<3>());
                    valid_systems.push_back(kv.first);
                    valid_sigma0.push_back(r.sigma0);
                }
            }

            // 系统间一致性检查
            if (valid_solutions.size() > 1) {
                Vector3d mean = Vector3d::Zero();
                for (const auto& sol : valid_solutions) {
                    mean += sol;
                }
                mean /= valid_solutions.size();
                
                double max_diff = 0;
                for (const auto& sol : valid_solutions) {
                    double diff = (sol - mean).norm();
                    max_diff = max(max_diff, diff);
                }
                
                if (max_diff > 10.0) {
                    cout << "历元 " << time_info.epoch_num << ": 系统间差异过大 (" << max_diff << "m)" << endl;
                    
                    // 找出最佳系统（sigma0最小）
                    int best_index = 0;
                    for (int i = 1; i < valid_sigma0.size(); i++) {
                        if (valid_sigma0[i] < valid_sigma0[best_index]) {
                            best_index = i;
                        }
                    }
                    
                    // 只保留最佳系统
                    Vector3d best_solution = valid_solutions[best_index];
                    char best_sys = valid_systems[best_index];
                    
                    // 重新设置有效解
                    valid_solutions.clear();
                    valid_solutions.push_back(best_solution);
                    valid_systems.clear();
                    valid_systems.push_back(best_sys);
                    
                    cout << "  使用最佳系统: " << best_sys << " (sigma0=" << valid_sigma0[best_index] << "m)" << endl;
                }
            }

            // ★ 第二步：输出结果并计算平均坐标 ★
            Vector3d sumXYZ(0,0,0);
            int valid_sys_count = 0;

            // 先输出所有系统的结果
            for (auto &kv : sys_results) {
                char sysc = kv.first;
                SolveResult &r = kv.second;
                
                outfile << left << setw(6) << string(1, sysc)
                        << setw(15) << r.Parameter(0) << setw(15) << r.Parameter(1) 
                        << setw(15) << r.Parameter(2) << setw(12) << r.Parameter(3)
                        << setw(12) << r.sigma0 << "\n";
            }

            // 只使用经过一致性检查的有效解来计算平均
            for (int i = 0; i < valid_solutions.size(); i++) {
                if (isReasonableCoordinate(valid_solutions[i](0), valid_solutions[i](1), valid_solutions[i](2))) {
                    sumXYZ += valid_solutions[i];
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