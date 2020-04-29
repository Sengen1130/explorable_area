#include "max_explorable_are.h"

int main()
{
    Explorable explorable;

    std::ofstream fout("max_explorable_area.dat");
    if (!fout)
    {
        std::cout << "ファイルをオープンできませんでした" << std::endl;
    }

    for (double j = 0; j < explorable.max_j; j++)
    {
        for (double i = 0; i < explorable.max_i; i++)
        {
            explorable.L = sqrt(pow(explorable.H, 2) + pow(i * 0.1, 2) + pow(explorable.d - j * 0.1, 2));
            explorable.phi = asin(explorable.H / explorable.L);
            explorable.theta = asin(i * 0.1 / sqrt(pow(explorable.L, 2) - pow(explorable.H, 2)));
            explorable.Fx = explorable.F * cos(explorable.phi) * sin(explorable.theta);
            explorable.Fy = explorable.F * cos(explorable.phi) * cos(explorable.theta);
            explorable.Fz = explorable.F * sin(explorable.phi);
            explorable.N_l = (explorable.m * explorable.g - explorable.Fz) / 2 + explorable.r / explorable.l * explorable.Fx;
            explorable.N_r = (explorable.m * explorable.g - explorable.Fz) / 2 - explorable.r / explorable.l * explorable.Fx;
            explorable.f_lx = explorable.Fx / 2 + pow(explorable.Fx, 2) * explorable.r / ((explorable.m * explorable.g - explorable.Fz) * explorable.l);
            explorable.f_rx = explorable.Fx / 2 - pow(explorable.Fx, 2) * explorable.r / ((explorable.m * explorable.g - explorable.Fz) * explorable.l);
            explorable.tau_brmax = explorable.r * sqrt(pow(explorable.myu * explorable.N_r, 2) - pow(explorable.f_rx, 2));

            for (int k = 0; k < 100; k++)
            {
                explorable.tau_br[k] = explorable.tau_brmax * k / 99;
            }

            for (int h = 0; h < 100; h++)
            {
                if (explorable.tau_br[h] <= explorable.r * explorable.Fy)
                {
                    explorable.tau_brmax2 = explorable.tau_br[h];
                }
            }

            explorable.nb = (explorable.l / (2 * explorable.r)) * explorable.tau_brmax2;
            explorable.Nt = -explorable.F * explorable.t * cos(explorable.phi) * sin(explorable.theta); //車体傾きなし
            explorable.Nb = explorable.nb;
            explorable.Nr = explorable.Nt + explorable.Nb;

            if (explorable.Nr > 0)
            {
                fout << i * 0.1 << " " << j * 0.1 << " " << 1 << std::endl;
            }
            else
            {
                fout << i * 0.1 << " " << j * 0.1 << " " << 0 << std::endl;
            }
        }
    }

    return 0;
}