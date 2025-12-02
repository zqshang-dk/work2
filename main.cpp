#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

// 单颗卫星的观测数据结构
struct SatObs {
    string prn;     // 卫星编号，如 C01
    double X;       // 卫星坐标 X (m)
    double Y;       // 卫星坐标 Y (m)
    double Z;       // 卫星坐标 Z (m)
    double P;       // 伪距 (m)
    double sigma;   // 观测噪声标准差 σ
};

// 单个历元的数据结构
struct EpochData {
    int epoch;      // 历元编号
    int gpsWeek;    // GPS周
    double gpsSec;  // 周内秒
    int satNum;     // 当前历元卫星数
    vector<SatObs> sats; // 卫星观测数据列表
};

// ----------------------------------------------------------------------
// 读取数据文件，解析所有历元
// ----------------------------------------------------------------------
vector<EpochData> readObsFile(const string& filename)
{
    vector<EpochData> allEpochs;
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "无法打开文件：" << filename << endl;
        return allEpochs;
    }

    string line;

    while (getline(fin, line)) {
        if (line.size() < 2) continue;

        if (line[0] == '#') {
            // 解析历元头 "#  1 2170 172800.000 25"
            EpochData epoch;
            stringstream ss(line.substr(1));

            ss >> epoch.epoch >> epoch.gpsWeek >> epoch.gpsSec >> epoch.satNum;

            epoch.sats.clear();

            // 读取该历元所有卫星记录
            for (int i = 0; i < epoch.satNum; i++) {
                getline(fin, line);
                stringstream sl(line);

                SatObs s;
                sl >> s.prn >> s.X >> s.Y >> s.Z >> s.P >> s.sigma;

                epoch.sats.push_back(s);
            }

            allEpochs.push_back(epoch);
        }
    }

    return allEpochs;
}


void LS(double X,double Y,double Z,double dT){
    double x = 0, y = 0, z = 0, dT = 0;
    double rou0 = sqrt((X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z));

    
}
// ----------------------------------------------------------------------
// 示例主函数：只是读取 + 打印前几行
// ----------------------------------------------------------------------
int main()
{
    string filename = "E:/STUDY/Sophomore1/最优估计/第二次上机实习/work2/CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt";

    auto epochs = readObsFile(filename);

    cout << "共读取历元数：" << epochs.size() << endl;

    // 打印第一历元前6颗卫星
    if (!epochs.empty()) {
        const auto& ep = epochs[0];

        cout << "第 " << ep.epoch << " 历元, 卫星数：" << ep.satNum << endl;

        for (int i = 0; i < min(6, ep.satNum); i++) {
            const auto& s = ep.sats[i];
            cout << s.prn << " "
                 << s.X << " " << s.Y << " " << s.Z << " "
                 << s.P << "  sigma=" << s.sigma << endl;
        }
    }

    return 0;
}
