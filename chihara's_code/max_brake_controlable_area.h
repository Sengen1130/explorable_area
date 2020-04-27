#include <iostream>
#include <fstream>
#include <cmath>
#define PI 3.14
using namespace std;

class Brake
{
public:
    const double m = 4.2;    //ローバ質量
    const double g = 9.8;    //重力加速度
    const double l = 0.286;  //トレッド間距離
    const double r = 0.092;  //タイヤ半径
    const double t = 0.063;  //タイヤシャフト中心とテザー位置間距離
    const double myu = 0.76; //動摩擦力
    double F;                //牽引張力
    double tether_xy_degree; //水平面内のテザー角度
    double tether_xz_degree; //水平面と垂直方向のテザー角度
    double phi;
    double theta;
    double Fx;
    double Fy;
    double Fz;
    double N_l;
    double N_r;
    double f_lx;
    double f_rx;
    double Nt;
    double Nb;
    double Nr;
    double tau_blmax; //左タイヤの最大ブレーキ力
    double tau_brmax; //右タイヤの最大ブレーキ力

    void set_parameter();
    void check_calulate_force();
    double degree_to_radian(double degree);
};