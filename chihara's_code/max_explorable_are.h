#include <iostream>
#include <fstream>
#include <cmath>
#define PI 3.14

class Explorable
{
public:
    const double m = 4.2;    //ローバーの質量[kg]
    const double g = 9.8;    //重力加速度[kg/s^2]
    const double l = 0.286;  //トレッド間距離
    const double r = 0.092;  //タイヤ半径
    const double t = 0.063;  //タイヤシャフト中心とテザー位置間距離
    const double myu = 0.76; //動摩擦係数
    double F = 20;           //リーダからの牽引力
    double d = 4;            //ウインチ間距離[m]
    double H = 1.5;          //ポール高さ[m]
    double max_i = d / 0.1;  //分割間隔
    double max_j = d / 0.1;  //分割間隔
    double theta;            //水平面内のテザー角度
    double phi;              //水平面と垂直方向のテザー角度
    double Fx;
    double Fy;
    double Fz;
    double N_l;
    double N_r;
    double f_lx;
    double f_rx;
    double tau_blmax; //左タイヤブレーキトルク最大値
    double tau_brmax; //右タイヤブレーキトルク最大値
    double tau_brmax2;
    double xF;
    double yF;
    double L; //ウインチのワイヤー長さ
    double Nt;
    double Nb;
    double Nr;
    double nb;
    double tau_br[100] = {};
};