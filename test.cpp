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
    Vector4d Parameter;
    double sigma0;  //验后单位权中误差
    Matrix4d cov;   //验后参数协方差矩阵
    bool success;   //是否解算成功
};

SolveResult pntpos(const vector<SatData>& obs, const EpochTime& time) {
    SolveResult res;
    res.success = false;

    vector<SatData> cleaned_obs;
    for (int i = 0; i < obs.size(); i++) {
        if (obs[i].pseudorange > 10000.0) {  // 只保留伪距大于10000米的卫星
            cleaned_obs.push_back(obs[i]);
        }
    }
    //太悲伤了！！！！！因为没有进行数据清洗导致02和03的粗差没有消除！
    
    int n = cleaned_obs.size();
    
    // 如果卫星数少于4颗，无法定位
    if (n < 4) {
        return res;
    }

    //数据初始化
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
        //用Eigen来定义矩阵和向量
        
        
        //填充H,l,W，此时需要遍历所有星历才可以实现
        for (int i = 0; i < n;i++){
            double P_obs = cleaned_obs[i].pseudorange;
            double Xs = cleaned_obs[i].sat_x;
            double Ys = cleaned_obs[i].sat_y;
            double Zs = cleaned_obs[i].sat_z; 
            double Var = cleaned_obs[i].variance;

            if(Var<=0)
                Var = 1.0;
            
            double P0 = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));
            
            H(i, 0) = (Xr - Xs) / P0;
            H(i, 1) = (Yr - Ys) / P0;
            H(i, 2) = (Zr - Zs) / P0;
            H(i, 3) = 1.0;

            //OMC残差
            l(i) = P_obs - (P0 + dt);

            //构造权矩阵
            W(i, i) = 1.0 / Var;
        }
        VectorXd dx = (H.transpose() * W * H).ldlt() .solve (H.transpose() * W * l);

        Xr += dx(0);
        Yr += dx(1);
        Zr += dx(2);
        dt += dx(3);

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

        double dist = sqrt((Xs - Xr) * (Xs - Xr) + (Ys - Yr) * (Ys - Yr) + (Zs - Zr) * (Zs - Zr));

        double P_cal = dist + dt;
        V(i) = P_cal - P_obs;
    }
    //计算验后单位权方差
    double vtwv = V.transpose() * W * V;
    double sigma0 = sqrt(vtwv / (n - 4));

    //计算验后协方差矩阵
    MatrixXd Qxx= (H.transpose() * W * H).inverse();

    Matrix4d Dxx = Qxx * sigma0 * sigma0;

    res.Parameter << Xr, Yr, Zr, dt;
    res.sigma0 = sigma0;
    res.cov = Dxx;
    res.success = true;

    return res;
}

// ECEF → ENU 转换（输入坐标为 m）
Vector3d ecef2enu(double X, double Y, double Z,
                  double X0, double Y0, double Z0)
{
    // 计算参考点的经纬度
    double dx = X - X0;
    double dy = Y - Y0;
    double dz = Z - Z0;

    double lon = atan2(Y0, X0);
    double p = sqrt(X0*X0 + Y0*Y0);
    double lat = atan2(Z0, p);

    // 构建旋转矩阵
    Matrix3d R;
    R << -sin(lon),              cos(lon),               0,
         -sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat),
          cos(lat)*cos(lon),  cos(lat)*sin(lon),  sin(lat);

    Vector3d enu = R * Vector3d(dx, dy, dz);
    return enu;
}

struct EvalResult {
    double meanX, meanY, meanZ;
    double stdX, stdY, stdZ;
    double meanE, meanN, meanU;
    double rmsE, rmsN, rmsU;
};

EvalResult evaluate(const vector<Vector4d>& allPos,
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

    // ================================
    // 3) 计算每个历元的 ENU
    // ================================
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
    double sumE = 0, sumN = 0, sumU = 0;
    for (int i = 0; i < m; i++) {
        sumE += E[i];
        sumN += N[i];
        sumU += U[i];
    }
    R.meanE = sumE / m;
    R.meanN = sumN / m;
    R.meanU = sumU / m;
    
    // 计算RMS
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
    vector<Vector4d> all_results;  // 存储每个历元的 X Y Z dt

    // 原始文件路径 (使用了 R"()" 语法，不需要双斜杠)
    string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt)";
    string outfile_path = "outpos_total.txt";

    ifstream infile(filename);
    ofstream outfile(outfile_path);

    if (!infile.is_open()) {
        cerr << "错误：无法打开文件 " << filename << endl;
        return -1;
    }

    // 设置输出精度为小数点后3位
    outfile << fixed << setprecision(3);

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
                
                SolveResult res = pntpos(epoch_obs, time_info);

                if (res.success) {
                    // ==========================================
                    //  按照图片格式输出
                    // ==========================================
                    
                    // 第一行：历元编号 GPS周 周秒 观测值数
                    outfile << "# " << left << setw(5) << time_info.epoch_num 
                            << setw(8) << time_info.week 
                            << setw(12) << fixed << setprecision(4) << time_info.tow 
                            << num_sats << endl;

                    // 第二行：表头
                    outfile << left << setw(15) << "X (m)" 
                            << setw(15) << "Y (m)" 
                            << setw(15) << "Z (m)" 
                            << "T(m)" << endl;

                    // 第三行：解算数值
                    outfile << left << fixed << setprecision(4) 
                            << setw(15) << res.Parameter(0) 
                            << setw(15) << res.Parameter(1) 
                            << setw(15) << res.Parameter(2) 
                            << res.Parameter(3) << endl;

                    // 第四行：验后单位权中误差
                    outfile << "验后单位权中误差: " << res.sigma0 << " (m)" << endl;

                    // 第五行：验后估计方差标题
                    outfile << "验后估计方差(m^2)" << endl;

                    // 第六至九行：4x4 协方差矩阵
                    for (int r = 0; r < 4; ++r) {
                        outfile << left << fixed << setprecision(4)
                                << setw(15) << res.cov(r, 0)
                                << setw(15) << res.cov(r, 1)
                                << setw(15) << res.cov(r, 2)
                                << setw(15) << res.cov(r, 3) << endl;
                    }

                    // 分隔线
                    outfile << "------------------------------------------------------------" << endl;

                    // 保存历元结果
                    all_results.push_back(res.Parameter);
                }
            }
        }
    }


    infile.close();
    outfile.close();
    cout << "处理完成，结果已保存至 " << outfile_path << endl;

    // 如果收集到了结果，进行精度评估
    if (!all_results.empty()) {
        // 设置参考坐标（可以设为第一个历元的解或真实坐标）
        double X0 = -1132914.6;  // 从你的第一个历元结果
        double Y0 = 6092528.9;
        double Z0 = 1504633.0;
        
        EvalResult eval = evaluate(all_results, X0, Y0, Z0, "result");
        
        // 输出评估结果
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