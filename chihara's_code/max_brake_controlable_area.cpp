#include "max_brake_controlable_area.h"
using namespace std;

void Brake::set_parameter()
{
    cout << "F=";
    cin >> F;
    cout << "tether_xy_degree=";
    cin >> tether_xy_degree;
    cout << "tether_xz_degree=";
    cin >> tether_xz_degree;
}

void Brake::check_calulate_force()
{
    cout << "N_l=" << N_l << endl;
    cout << "N_r=" << N_r << endl;
    cout << "f_lx=" << f_lx << endl;
    cout << "f_rx=" << f_rx << endl;
    cout << "tau_blmax=" << tau_blmax << endl;
    cout << "tau_brmax=" << tau_brmax << endl;
    cout << "phi=" << phi / PI * 180 << endl;
}

double Brake::degree_to_radian(double degree)
{
    return degree * PI / 180;
}

int main()
{
    Brake brake;

    ofstream fout("max_brake_controlable_area.dat");
    if (!fout)
    {
        cout << "ファイルをオープンできませんでした" << endl;
    }

    brake.set_parameter();

    brake.theta = brake.degree_to_radian(brake.tether_xy_degree); //水平面のテザー角度
    brake.phi = brake.degree_to_radian(brake.tether_xz_degree);   //水平面と垂直方向のテザー角度

    double tau_bl[500] = {};
    double tau_br[500] = {};
    double fb[500] = {};
    double nb[500] = {};

    brake.Fx = brake.F * cos(brake.phi) * sin(brake.theta);
    brake.Fy = brake.F * cos(brake.phi) * cos(brake.theta);
    brake.Fz = brake.F * sin(brake.phi);
    brake.N_l = (brake.m * brake.g - brake.Fz) / 2 + brake.r / brake.l * brake.Fx;
    brake.N_r = (brake.m * brake.g - brake.Fz) / 2 - brake.r / brake.l * brake.Fx;
    brake.f_lx = brake.Fx / 2 + pow(brake.Fx, 2) * brake.r / ((brake.m * brake.g - brake.Fz) * brake.l);
    brake.f_rx = brake.Fx / 2 - pow(brake.Fx, 2) * brake.r / ((brake.m * brake.g - brake.Fz) * brake.l);
    brake.tau_blmax = brake.r * sqrt(pow(brake.myu * brake.N_l, 2) - pow(brake.f_lx, 2));
    brake.tau_brmax = brake.r * sqrt(pow(brake.myu * brake.N_r, 2) - pow(brake.f_rx, 2));
    brake.check_calulate_force();

    for (int i = 0; i < 500; i++)
    {
        tau_br[i] = brake.tau_brmax * i / 499;
        tau_bl[i] = brake.tau_blmax * i / 499;
    }

    double nb_max[500] = {};

    for (int i = 0; i < 500; i++)
    {
        for (int j = 0; j < 500; j++)
        {
            if (tau_br[i] + tau_bl[j] <= brake.r * brake.Fy)
            {
                fb[i] = (1 / brake.r) * tau_br[i] + (1 / brake.r) * tau_bl[j];
                nb[j] = (brake.l / (2 * brake.r)) * tau_br[i] - (brake.l / (2 * brake.r)) * tau_bl[j];
                fout << fb[i] << "  " << nb[j] << endl;

                /*decide nb_max at roop i*/
                if (brake.tether_xy_degree > 0)
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

    double nb_Max = 0;

    if (brake.tether_xy_degree > 0)
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

    brake.Nt = -brake.F * brake.t * cos(brake.phi) * sin(brake.theta);
    brake.Nb = nb_Max;
    brake.Nr = brake.Nt + brake.Nb;
    cout << "nt=" << brake.Nt << endl;
    cout << "nr=" << brake.Nr << endl;
    cout << "------------------------" << endl;

    return 0;
}
