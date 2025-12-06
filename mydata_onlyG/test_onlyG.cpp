// //说明：根据已有的代码我们发现，E（Galileo）卫星数量较少，而且中误差都一样（均为4.0）
// // 所以这里我们考虑去除Galileo卫星，再次计算，看看结果有多大的差距
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <vector>
// #include <sstream>
// #include <iomanip>
// #include <cmath>
// #include <Eigen/Dense>

// using namespace std;
// using namespace Eigen;

// //光速（m/s）
// const double C_LIGHT = 299792458.0;
// const double A = 6378137.0;
// const double F = 1.0 / 298.257223563;
// const double E2 = 2 * F - F * F;

// // 定义一个结构体用来存储单颗卫星的观测数据
// struct SatData {
//     string id;      // 卫星ID (例如 C01)
//     double sat_x;        // 卫星 X 坐标 (m)
//     double sat_y;        // 卫星 Y 坐标 (m)
//     double sat_z;        // 卫星 Z 坐标 (m)
//     double pseudorange;  // 伪距观测值 (m) (对应文件第5列)
//     double variance;     // 误差方差 (对应文件第6列)
// };

// // 定义一个结构体存储当前历元的时间信息
// struct EpochTime {
//     int epoch_num; // 历元编号
//     int week;      // GPS周
//     double tow;    // 周内秒 (Time of Week)
// };

// //定义解算结果结构体
// struct SolveResult{
//     VectorXd Parameter;
//     double sigma0;  //验后单位权中误差
//     //Matrix4d cov;   //验后参数协方差矩阵
//     bool success;   //是否解算成功
// };

// // ECEF → BLH 转换函数（高精度迭代法）
// void ecef2blh(double X, double Y, double Z, double& B, double& L, double& H) {
//     double p = sqrt(X*X + Y*Y);
//     L = atan2(Y, X);                    // 经度（rad）
//     double B0 = atan2(Z, p * (1 - E2)); // 初始纬度
//     double B_prev;
//     double N;
//     do {
//         B_prev = B0;
//         N = A / sqrt(1 - E2 * sin(B0)*sin(B0));
//         H = p / cos(B0) - N;
//         B0 = atan2(Z, p * (1 - E2 * N / (N + H)));
//     } while (fabs(B0 - B_prev) > 1e-12);

//     B = B0;
//     N = A / sqrt(1 - E2 * sin(B)*sin(B));
//     H = p / cos(B) - N;

//     B = B * 180.0 / M_PI;   // 转为度
//     L = L * 180.0 / M_PI;
// }

// SolveResult pntpos_multi(const vector<SatData>& obs) {
//     SolveResult res;
//     res.success = false;
//     res.Parameter = VectorXd::Zero(4);          //// X,Y,Z, dtG
//     vector<SatData> cleaned_obs;
//     for (int i = 0; i < obs.size(); i++) {
//         if (obs[i].pseudorange > 10000000.0 && obs[i].variance>0&& obs[i].id[0] != 'E') {  // 只保留伪距大于10000米的卫星
//             cleaned_obs.push_back(obs[i]);
//         }
//     }
//     //太悲伤了！！！！！因为没有进行数据清洗导致02和03的粗差没有消除！
    
//     int n = cleaned_obs.size();
    
//     // 如果卫星数少于6颗，无法定位
//     if (n < 4) {
//         return res;
//     }
   
//     //数据初始化（使用传入的初始值）
//     double Xr = 0.0;
//     double Yr = 0.0;
//     double Zr = 0.0;
//     double dt_GPS = 0.0;
//     //double dt_BDS = 0.0;
//     //double dt_Gal = 0.0;
//     //double dt_GLO = 0.0;

//     //迭代参数
//     const int maxIter = 10;
//     const double eps = 1e-4;  //收敛阈值（m）
    
//     MatrixXd H(n, 4);  //设计矩阵
//     VectorXd l(n);     //OMC观测值
//     MatrixXd W = MatrixXd::Zero(n, n);  //创建一个n*n的零矩阵

//     for (int iter = 0; iter < maxIter; ++iter){
//         //用Eigen来定义矩阵和向量    
//         //填充H,l,W，此时需要遍历所有星历才可以实现
//         for (int i = 0; i < n;i++){
//             double P_obs = cleaned_obs[i].pseudorange;
//             double Xs = cleaned_obs[i].sat_x;
//             double Ys = cleaned_obs[i].sat_y;
//             double Zs = cleaned_obs[i].sat_z; 
//             double Var = cleaned_obs[i].variance;
//             //double P0 = sqrt((Xs - x(0)) * (Xs - x(0)) + (Ys - x(1)) * (Ys - x(1)) + (Zs - x(2)) * (Zs - x(2)));
            
            
//             if(Var<=0)
//                 Var = 1.0;
            
//             double P0 = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
//             //H.row(i).setZero();
            
            
//             H(i, 0) = (Xr - Xs) / P0;
//             H(i, 1) = (Yr - Ys) / P0;
//             H(i, 2) = (Zr - Zs) / P0;
//             // H(i, 3) = 1.0;

//             //OMC残差
//             // l(i) = P_obs - (P0 + dt);

//             //钟差列：根据系统设置1
//             char sys = cleaned_obs[i].id[0];
//             if(sys=='G'){
//                 H(i, 3) = 1.0;
//             }
                
//             // else if(sys=='C'){
//             //     H(i, 3) = 0;
//             //     H(i, 4) = 1.0;
//             //     H(i, 5) = 0;

//             // }
//             // else if(sys=='E'){
//             //     H(i, 3) = 0;
//             //     H(i, 4) = 0;
//             //     H(i, 5) = 1.0;
//             //     H(i, 6) = 0;
//             // }
//             // else if(sys=='R'){
//             //     H(i, 3) = 0;
//             //     H(i, 4) = 0;
//             //     H(i, 5) = 1.0;

//             // }
//             // else
//             //     H(i, 3) = 1.0;  //如果没写的话默认是GPS卫星
            
            
//             double dt_sys = 0.0;
//             if (sys == 'G')
//                 dt_sys = dt_GPS;
//             // else if (sys == 'C')
//             //     dt_sys = dt_BDS;
//             // else if (sys == 'E')
//             //     dt_sys = dt_Gal;
//             // else if (sys == 'R')
//             //     dt_sys = dt_GLO;

//             //OMC残差
//             l(i) = P_obs - (P0 + dt_sys);

//             //构造权矩阵
//             W(i, i) = 1.0 / Var;
//         }
//         VectorXd dx = (H.transpose() * W * H).ldlt() .solve (H.transpose() * W * l);

//         Xr += dx(0);
//         Yr += dx(1);
//         Zr += dx(2);
//         dt_GPS+=dx(3);
//         //dt_BDS+=dx(4);
//         //dt_Gal+=dx(5);
//         //dt_GLO += dx(5);

//         double pos_delta = sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2));
//         if(pos_delta<eps){
//             break;
//         }
//     }

//     //精度评定
//     VectorXd V(n);  //残差
//     for (int i = 0; i < n;i++){
//         double P_obs = cleaned_obs[i].pseudorange;
//         double Xs = cleaned_obs[i].sat_x;
//         double Ys = cleaned_obs[i].sat_y;
//         double Zs = cleaned_obs[i].sat_z;
//         char sys = cleaned_obs[i].id[0];
        
//         double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
//         double dt_sys = dt_GPS;
//         double P_cal = dist + dt_sys;
//         V(i) = P_cal - P_obs;
//     }
//     //计算验后单位权方差
//     double vtwv = V.transpose() * W * V;
//     double sigma0 = sqrt(vtwv / (n - 4));

//     res.Parameter << Xr, Yr, Zr, dt_GPS ;
//     res.sigma0 = sigma0;

//     res.success = true;

//     return res;
// }



// int main() {
//     vector<VectorXd> all_results;  // 存储每个历元的 X Y Z dt

//     // 原始文件路径 (使用了 R"()" 语法，不需要双斜杠)
//     string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\mydata\2059329.25ObsCorr)";
//     string outfile_path = "outpos_total.txt";
//     string blhfile_path = "blh_results.txt";

//     ifstream infile(filename);
//     ofstream outfile(outfile_path);
//     ofstream blhfile(blhfile_path);  // 新增

//     if (!infile.is_open()) {
//         cerr << "错误：无法打开文件 " << filename << endl;
//         return -1;
//     }
     
//     if (!blhfile.is_open()) {
//         cerr << "错误：无法创建BLH输出文件" << endl;
//         return -1;
//     }

//     // 设置输出精度为小数点后3位
//     outfile << fixed << setprecision(3);
//     blhfile << fixed << setprecision(8); 

//     blhfile << "# Epoch Week TOW avg_lat(deg) avg_lon(deg) avg_h(m)\n";


//     string line;
//     while (getline(infile, line)) {
//         // 跳过空行
//         if (line.empty()) continue;

//         // 检查是否是历元头 (以 '#' 开头)
//         if (line[0] == '#') {
//             EpochTime time_info;
//             int num_sats = 0;
//             char hash_char; 

//             // 解析头文件: "#  1 2170 172800.000 25"
//             stringstream ss_header(line);
//             ss_header >> hash_char >> time_info.epoch_num >> time_info.week >> time_info.tow >> num_sats;

//             // 读取该历元下的所有卫星数据
//             vector<SatData> epoch_obs;
//             epoch_obs.reserve(num_sats); 

//             for (int i = 0; i < num_sats; ++i) {
//                 if (!getline(infile, line)) break;
                
//                 stringstream ss_obs(line);
//                 SatData sat;
//                 // 解析每行卫星数据
//                 ss_obs >> sat.id >> sat.sat_x >> sat.sat_y >> sat.sat_z >> sat.pseudorange >> sat.variance;
                
//                 epoch_obs.push_back(sat);
//             }

//             // 进行定位解算
//             if (epoch_obs.size() >= 6) {
//                 vector<double> result(4, 0.0); // 存放解 [x, y, z, dt(4)]
                
//                 //SolveResult res = pntpos_multi(epoch_obs);

//                 SolveResult res=pntpos_multi(epoch_obs);
                
//                 vector<VectorXd> all_positions;

//                 if (res.success) {
//                     //all_results.push_back(res.Parameter);

//                     double X = res.Parameter(0);
//                     double Y = res.Parameter(1);
//                     double Z = res.Parameter(2);
//                     double Lat, Lon, H;
//                     ecef2blh(X, Y, Z, Lat, Lon, H);

                   
//                     // ===== 输出头部 =====
//                     outfile << "# " << setw(4) << time_info.epoch_num
//                             << setw(8) << time_info.week
//                             << setw(14) << fixed << setprecision(4) << time_info.tow
//                             << setw(4) << epoch_obs.size() << "\n";

//                     // ===== 输出标题行 =====
//                     outfile << "X (m)          Y (m)          Z (m)          "
//                             << "dT_GPS(m)      \n";

//                     // ===== 输出数据行 =====
//                     outfile << fixed << setprecision(4)
//                             << setw(15) << X
//                             << setw(15) << Y
//                             << setw(15) << Z
//                             << setw(15) << res.Parameter(3)   // dT_GPS
//                             //<< setw(15) << res.Parameter(4)   // dT_BDS
//                             //<< setw(15) << res.Parameter(5)   // dT_GAL
//                             //<< setw(15) << res.Parameter(5)   // dT_GLO
//                             << "\n\n";
//                     outfile << " sigma0: " << fixed << setprecision(3) << res.sigma0 << " m\n\n";
//                     outfile << "------------------------------------------------------------" << endl;
                    
//                     // ===== 新增：输出到独立的BLH文件 =====
//                     blhfile << setw(5) << time_info.epoch_num  // Epoch
//                             << setw(6) << time_info.week       // Week
//                             << setw(14) << fixed << setprecision(4) << time_info.tow  // TOW
//                             << "  "  // 两个空格分隔
//                             << setw(15) << fixed << setprecision(8) << Lat    // Latitude
//                             << " "
//                             << setw(15) << fixed << setprecision(8) << Lon    // Longitude
//                             << " "
//                             << setw(12) << fixed << setprecision(4) << H      // Height
//                             << "\n";

//                     // 保存历元结果
//                     all_results.push_back(res.Parameter);
//                 }
//             }
//         }
//     }


//     infile.close();
//     outfile.close();
//     blhfile.close();
//     cout << "处理完成，结果已保存至 " << outfile_path << endl;

//     return 0;
// }

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
    res.Parameter = VectorXd::Zero(4);          //// X,Y,Z, dtG, dtC, dtE, dtR

    vector<SatData> cleaned_obs_initial;
    for (int i = 0; i < obs.size(); i++) {
        if (obs[i].pseudorange > 10000000.0 && obs[i].variance>0) {  // 只保留伪距大于10000米的卫星
            cleaned_obs_initial.push_back(obs[i]);
        }
    }

     // 如果卫星数少于4颗，无法定位
    if (cleaned_obs_initial.size() < 4) {
        return res;
    }
    //太悲伤了！！！！！因为没有进行数据清洗导致02和03的粗差没有消除！
    

    //又悲伤了！没有用3-sigma剔除异常值，算出来的还有粗差点！
    //vector<SatData> cleaned_obs = cleaned_obs_initial;


    //int n = cleaned_obs_initial.size();

    //---------------------------------------------------------
    // 2) ========== 三倍中误差（3σ）粗差剔除：一次性处理 ========== 
    //---------------------------------------------------------
    // 先用初值计算一次残差
    double X0 = 0, Y0 = 0, Z0 = 0;
    vector<double> residuals;
    residuals.reserve(cleaned_obs_initial.size());

    for (auto& s : cleaned_obs_initial) {
        double P0 = sqrt((s.sat_x - X0)*(s.sat_x - X0) +
                         (s.sat_y - Y0)*(s.sat_y - Y0) +
                         (s.sat_z - Z0)*(s.sat_z - Z0));
        residuals.push_back(s.pseudorange - P0);
    }

    // 计算残差均值与标准差
    double mean = 0, sigma = 0;
    for (double r : residuals) mean += r;
    mean /= residuals.size();

    for (double r : residuals) sigma += (r - mean) * (r - mean);
    sigma = sqrt(sigma / residuals.size());

    // 再根据 3σ 剔除粗差
    vector<SatData> cleaned_obs;
    for (int i = 0; i < cleaned_obs_initial.size(); i++) {
        if (fabs(residuals[i] - mean) <= 3 * sigma)  // <-- 3σ判断
            cleaned_obs.push_back(cleaned_obs_initial[i]);
    }
   
    int n = cleaned_obs.size();
    if (n < 4) return res;
    
    //数据初始化（使用传入的初始值）
    double Xr = 0.0;
    double Yr = 0.0;
    double Zr = 0.0;
    double dt_GPS = 0.0;
    // double dt_BDS = 0.0;
    // double dt_Gal = 0.0;
    // double dt_GLO = 0.0;

    //迭代参数
    const int maxIter = 10;
    const double eps = 1e-4;  //收敛阈值（m）
    
    MatrixXd H(n, 4);  //设计矩阵
    VectorXd l(n);     //OMC观测值
    MatrixXd W = MatrixXd::Zero(n, n);  //创建一个n*n的零矩阵

    for (int iter = 0; iter < maxIter; ++iter){
        //用Eigen来定义矩阵和向量    
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
            
            
            H(i, 0) = (Xr - Xs) / P0;
            H(i, 1) = (Yr - Ys) / P0;
            H(i, 2) = (Zr - Zs) / P0;
            // H(i, 3) = 1.0;

            //OMC残差
            // l(i) = P_obs - (P0 + dt);

            //钟差列：根据系统设置1
            char sys = cleaned_obs[i].id[0];
            if(sys=='G'){
                H(i, 3) = 1.0;
            }
                
            // else if(sys=='C'){
            //     H(i, 3) = 0;
            //     H(i, 4) = 1.0;
            //     H(i, 5) = 0;
            //     H(i, 6) = 0;
            // }
            // else if(sys=='E'){
            //     H(i, 3) = 0;
            //     H(i, 4) = 0;
            //     H(i, 5) = 1.0;
            //     H(i, 6) = 0;
            // }
            // else if(sys=='R'){
            //     H(i, 3) = 0;
            //     H(i, 4) = 0;
            //     H(i, 5) = 0;
            //     H(i, 6) = 1.0;
            // }
            // else
            //     H(i, 3) = 1.0;  //如果没写的话默认是GPS卫星
            
            
            double dt_sys = 0.0;
            if (sys == 'G')
                dt_sys = dt_GPS;
            // else if (sys == 'C')
            //     dt_sys = dt_BDS;
            // else if (sys == 'E')
            //     dt_sys = dt_Gal;
            // else if (sys == 'R')
            //     dt_sys = dt_GLO;

            //OMC残差
            l(i) = P_obs - (P0 + dt_sys);

            //构造权矩阵
            W(i, i) = 1.0 / Var;
        }
        VectorXd dx = (H.transpose() * W * H).ldlt() .solve (H.transpose() * W * l);

        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt_GPS+=dx(3);
        // dt_BDS+=dx(4);
        // dt_Gal+=dx(5);
        // dt_GLO += dx(6);

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
        
        double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
        double dt_sys = dt_GPS ;
        double P_cal = dist + dt_sys;
        V(i) = P_cal - P_obs;
    }
    //计算验后单位权方差
    double vtwv = V.transpose() * W * V;
    double sigma0 = sqrt(vtwv / (n - 4));

    res.Parameter << Xr, Yr, Zr, dt_GPS ;
    res.sigma0 = sigma0;

    res.success = true;

    return res;
}



int main() {
    vector<VectorXd> all_results;  // 存储每个历元的 X Y Z dt

    // 原始文件路径 (使用了 R"()" 语法，不需要双斜杠)
    string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\mydata\2059329.25ObsCorr)";
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
                vector<double> result(4, 0.0); // 存放解 [x, y, z, dt(4)]
                
                //SolveResult res = pntpos_multi(epoch_obs);

                SolveResult res=pntpos_multi(epoch_obs);
                
                vector<VectorXd> all_positions;

                if (res.success) {
                    //all_results.push_back(res.Parameter);

                    double X = res.Parameter(0);
                    double Y = res.Parameter(1);
                    double Z = res.Parameter(2);
                    double Lat, Lon, H;
                    ecef2blh(X, Y, Z, Lat, Lon, H);

                   
                    // ===== 输出头部 =====
                    outfile << "# " << setw(4) << time_info.epoch_num
                            << setw(8) << time_info.week
                            << setw(14) << fixed << setprecision(4) << time_info.tow
                            << setw(4) << epoch_obs.size() << "\n";

                    // ===== 输出标题行 =====
                    outfile << "X (m)          Y (m)          Z (m)          "
                            << "dT_GPS(m)      \n";

                    // ===== 输出数据行 =====
                    outfile << fixed << setprecision(4)
                            << setw(15) << X
                            << setw(15) << Y
                            << setw(15) << Z
                            << setw(15) << res.Parameter(3)   // dT_GPS
                            //<< setw(15) << res.Parameter(4)   // dT_BDS
                            //<< setw(15) << res.Parameter(5)   // dT_GAL
                            //<< setw(15) << res.Parameter(6)   // dT_GLO
                            << "\n\n";
                    outfile << " sigma0: " << fixed << setprecision(3) << res.sigma0 << " m\n\n";
                    outfile << "------------------------------------------------------------" << endl;
                    
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

    return 0;
}