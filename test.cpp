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

// 你的定位解算函数接口
// 参数：
//   obs: 当前历元所有可见卫星的数据列表
//   time: 当前历元的时间信息
//   out_pos: 输出参数，存放在此处 [x, y, z, dt]
void pntpos(const vector<SatData>& obs, const EpochTime& time, vector<double>& out_pos) {
    int n = obs.size();
    
    // 如果卫星数少于4颗，无法定位
    if (n < 4) {
        return;
    }

    // ==========================================
    //在此处编写最小二乘解算 (Least Squares)
    // ==========================================
    
    // 提示：如何获取第 i 颗卫星的数据 (假设循环变量是 i)
    // double P_obs = obs[i].pseudorange;  // 伪距观测值
    // double Xs    = obs[i].sat_x;        // 卫星X坐标
    // double Ys    = obs[i].sat_y;        // 卫星Y坐标
    // double Zs    = obs[i].sat_z;        // 卫星Z坐标
    // double Var   = obs[i].variance;     // 方差

    //数据初始化
    double Xr = 0;
    double Yr = 0;
    double Zr = 0;
    double dt = 0;

    //迭代参数
    const int maxIter = 10;
    const double eps = 1e-4;  //收敛阈值（m）

    for (int iter = 0; iter < maxIter; ++iter)
    {
        // double H[n][4];
        //     H(i, 0) = (Xr - Xs) / P0;
        //     H(i, 1) = (Yr - Ys) / P0;
        //     H(i, 2) = (Zr - Zs) / P0;
        //     H(i, 3) = 1;

        //用Eigen来定义矩阵和向量
        MatrixXd H(n, 4);  //设计矩阵n*4
        VectorXd l(n);     //OMC观测值
        MatrixXd W = MatrixXd::Zero(n, n);  //创建一个n*n的零矩阵

        //填充H,l,W，此时需要遍历所有星历才可以实现
        for (int i = 0; i < n;i++){
            double P_obs = obs[i].pseudorange;
            double Xs = obs[i].sat_x;
            double Ys = obs[i].sat_y;
            double Zs = obs[i].sat_z; 
            double Var = obs[i].variance;

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


        // for (int i = 0; i < n; i++)
        // {
            

            
        //     // 设计矩阵
        //     //  double H[25][4] = [((Xr - Xs) / P0, (Yr - Ys) / P0, (Zr - Zs) / P0, 1),];
            

        //     // OMC观测值计算
        //     double l[i] = P_obs - (P0 + c * dt);

        //     // 构造权矩阵
        //     double W[n][n];
        //     W[i][i] = 1.0 / var;

        //     double dx[4] = (H.transpose() * W * H).inverse() * (H.tranpose() * W * l[i]);
        //     Xr += dx(0);
        //     Yr += dx(1);
        //     Zr += dx(2);
        //     dt += dx(3);

        //     if (dx.norm() < 1e-4)
        //     {
        //         break;
        //     }
        // }
        

    // 模拟输出结果 (这只是为了演示，你需要算出真实值)
    out_pos.resize(4);
    out_pos[0] = Xr; // 用户 X
    out_pos[1] = Yr; // 用户 Y
    out_pos[2] = Zr; // 用户 Z
    out_pos[3] = dt; // 钟差
}

int main() {
    // 原始文件路径 (使用了 R"()" 语法，不需要双斜杠)
    string filename = R"(E:\STUDY\Sophomore1\最优估计\第二次上机实习\work2\CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt)";
    string outfile_path = "outpos_cpp.txt";

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
                
                pntpos(epoch_obs, time_info, result);

                // 输出格式: GPSWeek  GPSSec  UserX  UserY  UserZ
                outfile << time_info.week << " " 
                        << time_info.tow << " " 
                        << setw(14) << result[0] << " " 
                        << setw(14) << result[1] << " " 
                        << setw(14) << result[2] << endl;
            }
        }
    }

    infile.close();
    outfile.close();
    cout << "处理完成，结果已保存至 " << outfile_path << endl;

    return 0;
}