#include "max_brake_controlable_area.h"

void Brake::set_parameter()
{
    std::cout << "F=";
    std::cin >> F;
    std::cout << "tether_xy_degree=";
    std::cin >> tether_xy_degree;
    std::cout << "tether_xz_degree=";
    std::cin >> tether_xz_degree;
}

void Brake::check_calulate_force()
{
    std::cout << "N_l=" << N_l << std::endl;
    std::cout << "N_r=" << N_r << std::endl;
    std::cout << "f_lx=" << f_lx << std::endl;
    std::cout << "f_rx=" << f_rx << std::endl;
    std::cout << "tau_blmax=" << tau_blmax << std::endl;
    std::cout << "tau_brmax=" << tau_brmax << std::endl;
    std::cout << "phi=" << phi / PI * 180 << std::endl;
    std::cout << "theta=" << theta / PI * 180 << std::endl;
    std::cout << "nt=" << Nt << std::endl;
    std::cout << "nr=" << Nr << std::endl;
    std::cout << "nb_max=" << nb_Max << std::endl;
}

double Brake::degree_to_radian(double degree)
{
    return degree * PI / 180;
}

int main()
{
    Brake brake;

    std::ofstream fout("max_brake_controlable_area.dat");
    if (!fout)
    {
        std::cout << "ファイルをオープンできませんでした" << std::endl;
    }

    brake.set_parameter();

    brake.theta = brake.degree_to_radian(brake.tether_xy_degree); //水平面のテザー角度
    brake.phi = brake.degree_to_radian(brake.tether_xz_degree);   //水平面と垂直方向のテザー角度

    brake.Fx = brake.F * cos(brake.phi) * sin(brake.theta);
    brake.Fy = brake.F * cos(brake.phi) * cos(brake.theta);
    brake.Fz = brake.F * sin(brake.phi);
    brake.N_l = (brake.m * brake.g - brake.Fz) / 2 + brake.r / brake.l * brake.Fx;
    brake.N_r = (brake.m * brake.g - brake.Fz) / 2 - brake.r / brake.l * brake.Fx;
    brake.f_lx = brake.Fx / 2 + pow(brake.Fx, 2) * brake.r / ((brake.m * brake.g - brake.Fz) * brake.l);
    brake.f_rx = brake.Fx / 2 - pow(brake.Fx, 2) * brake.r / ((brake.m * brake.g - brake.Fz) * brake.l);
    brake.tau_blmax = brake.r * sqrt(pow(brake.myu * brake.N_l, 2) - pow(brake.f_lx, 2));
    brake.tau_brmax = brake.r * sqrt(pow(brake.myu * brake.N_r, 2) - pow(brake.f_rx, 2));

    for (int i = 0; i < 500; i++)
    {
        brake.tau_br[i] = brake.tau_brmax * i / 499;
        brake.tau_bl[i] = brake.tau_blmax * i / 499;
    }

    for (int i = 0; i < 500; i++)
    {
        for (int j = 0; j < 500; j++)
        {
            if (brake.tau_br[i] + brake.tau_bl[j] <= brake.r * brake.Fy)
            {
                brake.fb[i] = (1 / brake.r) * brake.tau_br[i] + (1 / brake.r) * brake.tau_bl[j];
                brake.nb[j] = (brake.l / (2 * brake.r)) * brake.tau_br[i] - (brake.l / (2 * brake.r)) * brake.tau_bl[j];
                fout << brake.fb[i] << "  " << brake.nb[j] << std::endl;

                /*decide nb_max at roop i*/
                if (brake.tether_xy_degree > 0)
                {
                    if (brake.nb[j] > brake.nb_max[i])
                        brake.nb_max[i] = brake.nb[j];
                }
                else
                {
                    if (brake.nb[j] < brake.nb_max[i])
                        brake.nb_max[i] = brake.nb[j];
                }
            }
        }
    }

    if (brake.tether_xy_degree > 0)
    {
        for (int i = 0; i < 500; i++)
        {
            if (brake.nb_max[i] > brake.nb_Max)
                brake.nb_Max = brake.nb_max[i];
        }
    }
    else
    {
        for (int i = 0; i < 500; i++)
        {
            if (brake.nb_max[i] < brake.nb_Max)
                brake.nb_Max = brake.nb_max[i];
        }
    }

    brake.Nt = -brake.F * brake.t * cos(brake.phi) * sin(brake.theta);
    brake.Nb = brake.nb_Max;
    brake.Nr = brake.Nt + brake.Nb;

    brake.check_calulate_force();

    return 0;
}
