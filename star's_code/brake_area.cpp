/*実現可能なロボットにかかる力・モーメント（上からの力追加）*/
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#define PI 3.14

double deg2rad(double a)
{
  double b;
  b = a * PI / 180;
  return b;
}

int main()
{
  ofstream fout("kadai3.dat");
  if (!fout)
  {
    cout << "ファイルをオープンできませんでした" << endl;
    return 1;
  }
  const double m = 4.2;   //ローバ質量
  const double g = 9.8;   //重力加速度
  const double l = 0.286; //トレッド間距離
  const double r = 0.092; //タイヤ半径
  const double t = 0.063; //タイヤシャフト中心とテザー位置間距離
  const double myu = 0.76;
  //const double n = 0.7;
  //const  double tether_xz = 41;

  double F; //牽引張力
  cout << "F=";
  cin >> F;
  double tether_xy; //水平面内のテザー角度
  cout << "tether_xy=";
  cin >> tether_xy;
  double tether_xz; //水平面と垂直方向のテザー角度
  cout << "tether_xz=";
  cin >> tether_xz;

  //double phi = asin(m*g*(1-n)/F);
  double theta = deg2rad(tether_xy); //水平面のテザー角度
  double phi = deg2rad(tether_xz);   //水平面と垂直方向のテザー角度

  double Fx, Fy, Fz, N_l, N_r, f_lx, f_rx;
  double tau_blmax; //左タイヤの最大ブレーキ力
  double tau_brmax; //右タイヤの最大ブレーキ力
  double Nt, Nb, Nr;

  double tau_bl[500] = {};
  double tau_br[500] = {};
  double fb[500] = {};
  double nb[500] = {};

  Fx = F * cos(phi) * sin(theta);
  Fy = F * cos(phi) * cos(theta);
  Fz = F * sin(phi);
  N_l = (m * g - Fz) / 2 + r / l * Fx;
  N_r = (m * g - Fz) / 2 - r / l * Fx;
  f_lx = Fx / 2 + pow(Fx, 2) * r / ((m * g - Fz) * l);
  f_rx = Fx / 2 - pow(Fx, 2) * r / ((m * g - Fz) * l);
  tau_blmax = r * sqrt(pow(myu * N_l, 2) - pow(f_lx, 2));
  tau_brmax = r * sqrt(pow(myu * N_r, 2) - pow(f_rx, 2));

  /*確認用*/
  cout << "N_l=" << N_l << endl;
  cout << "N_r=" << N_r << endl;
  cout << "f_lx=" << f_lx << endl;
  cout << "f_rx=" << f_rx << endl;
  cout << "tau_blmax=" << tau_blmax << endl;
  cout << "tau_brmax=" << tau_brmax << endl;
  //cout << "total=" << r*Fy << endl;
  cout << "phi=" << phi / PI * 180 << endl;
  //cout << "tow=" << r*ft*cos(phi) << endl;

  for (int i = 0; i < 500; i++)
  {
    tau_br[i] = tau_brmax * i / 499;
    tau_bl[i] = tau_blmax * i / 499;
  }

  double nb_max[500] = {};

  for (int i = 0; i < 500; i++)
  {
    for (int j = 0; j < 500; j++)
    {
      if (tau_br[i] + tau_bl[j] <= r * Fy)
      {
        fb[i] = (1 / r) * tau_br[i] + (1 / r) * tau_bl[j];
        nb[j] = (l / (2 * r)) * tau_br[i] - (l / (2 * r)) * tau_bl[j];
        fout << fb[i] << "  " << nb[j] << endl;

        /*decide nb_max at roop i*/
        if (tether_xy > 0)
        {
          if (nb[j] > nb_max[i])
            nb_max[i] = nb[j];
        }
        else
        {
          if (nb[j] < nb_max[i])
            nb_max[i] = nb[j];
        }
      }
    }
  }

  /*  double fb_max = 0;

  for(int i = 0;i<100;i++)
    {
      if(fb[i] > fb_max) fb_max = fb[i];
    }
    cout << "fb_max=" << fb_max << endl;*/

  double nb_Max = 0;

  if (tether_xy > 0)
  {
    for (int i = 0; i < 500; i++)
    {
      if (nb_max[i] > nb_Max)
        nb_Max = nb_max[i];
    }
  }
  else
  {
    for (int i = 0; i < 500; i++)
    {
      if (nb_max[i] < nb_Max)
        nb_Max = nb_max[i];
    }
  }

  cout << "nb_max=" << nb_Max << endl;

  Nt = -F * t * cos(phi) * sin(theta); //車体傾きなし
  //double a = 45;
  //double alpha = deg2rad(a);
  //Nt = -Fx*t*cos(alpha); //車体傾きあり
  Nb = nb_Max;
  Nr = Nt + Nb;
  cout << "nt=" << Nt << endl;
  cout << "nr=" << Nr << endl;
  cout << "------------------------" << endl;

  return 0;
}
