/*実現可能なロボットにかかる力・モーメント（上からの力追加）*/
#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;
#define PI 3.14

double deg2rad(double a)
{
  double b;
  b=a*PI/180;
  return b;
}

int main()
{
  ofstream fout("limit_area.dat");
  if(!fout)
    {
      cout << "ファイルをオープンできませんでした" << endl;
      return 1;
    }
  const double m = 4.2;		//ローバーの質量[kg]
  const double g = 9.8;		//重力加速度[kg/s^2]
  const double l = 0.286;	//トレッド間距離
  const double r = 0.092; 	//タイヤ半径
  const double t = 0.063;	//タイヤシャフト中心とテザー位置間距離
  const double myu = 0.76;	//動摩擦係数
  //const double n = 0.7;
  //const  double tether_xz = 41;
  //const double F = 20;

  double F = 20;		//リーダからの牽引力
  //double phi = asin(m*g*(1-n)/F);
  double theta;		//水平面内のテザー角度
  double phi;		//水平面と垂直方向のテザー角度
  double Fx,Fy,Fz,N_l,N_r,f_lx,f_rx;
  double tau_blmax;		//左タイヤブレーキトルク最大値
  double tau_brmax;		//右タイヤブレーキトルク最大値
  double tau_brmax2;
  double xF,yF;
  double L;		//ウインチのワイヤー長さ
  double d=4;		//ウインチ間距離[m]
  double H=1.5;		//ポール高さ[m]
  double Nt,Nb,Nr,nb;
  double tau_br[100] = {};

  double i_max = d/0.1;		//分割間隔
  double j_max = d/0.1;		//分割間隔
  double i,j;

  for(j=0;j<j_max;j++)
    {
      for(i=0;i<i_max;i++)
	{
	  L = sqrt(pow(H,2)+pow(i*0.1,2)+pow(d-j*0.1,2));
	  phi = asin(H/L);
	  theta = asin(i*0.1/sqrt(pow(L,2)-pow(H,2)));
	  Fx = F*cos(phi)*sin(theta);
	  Fy = F*cos(phi)*cos(theta);
	  Fz = F*sin(phi); 
	  N_l = (m*g-Fz)/2+r/l*Fx;
	  N_r = (m*g-Fz)/2-r/l*Fx;
	  f_lx = Fx/2+pow(Fx,2)*r/((m*g-Fz)*l);
	  f_rx = Fx/2-pow(Fx,2)*r/((m*g-Fz)*l);
	  //tau_blmax = r*sqrt(pow(myu*N_l,2)-pow(f_lx,2));
	  tau_brmax = r*sqrt(pow(myu*N_r,2)-pow(f_rx,2));
  
	  /*確認用*/
	  /*
	  cout << "N_l=" << N_l << endl;
	  cout << "N_r=" << N_r << endl;
	  cout << "f_lx=" << f_lx << endl;
	  cout << "f_rx=" << f_rx << endl;
	  cout << "tau_blmax=" << tau_blmax << endl;
	  cout << "tau_brmax=" << tau_brmax << endl;
	  //cout << "total=" << r*Fy << endl;
	  cout << "phi=" << phi/PI*180 << endl;
	  //cout << "tow=" << r*ft*cos(phi) << endl; 
	  */
	  
	  for(int k=0;k<100;k++)
	    {
	      tau_br[k] = tau_brmax*k/99;
	    }

	  for(int h=0;h<100;h++)
	    {
	      if(tau_br[h] <= r*Fy)
		{
		  tau_brmax2 =tau_br[h];
		}
	    }
	  nb = (l/(2*r))*tau_brmax2;
	  Nt = -F*t*cos(phi)*sin(theta);//車体傾きなし
	  Nb = nb;
	  Nr = Nt + Nb;
	  if(Nr>0)
	    {
	      fout << i*0.1 << " " << j*0.1 << " " << 1 << endl;
	    }
	  else
	    {
	      fout << i*0.1 << " " << j*0.1 << " " << 0 << endl;
	    }
	}
    }

    return 0;
}


