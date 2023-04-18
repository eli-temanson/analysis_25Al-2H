/*
 * Simple energy  loss integrators for forward and reverse cases. 
 *
 * Dec 2022, Gordon McCann
 */
#include "catima/gwm_integrators.h"

namespace catima {

    double integrate_energyloss(Projectile& proj, const Material& mat, const Config& c)
    {
        static double s_estep_max = 0.001;
        static int s_depth_max = 100;
        int depth = 0;
        double e_in = proj.T; // MeV/u
        double e_final = e_in;
        double x_step = 0.25*mat.thickness(); //g/cm^2
        double x_traversed = 0.0;
        double e_step = dedx(proj, mat, c)*x_step;
		double A_recip = 1.0/proj.A;

        while(true)
        {
            if(e_step/e_final > s_estep_max && depth < s_depth_max)
            {
                ++depth;
                x_step *= 0.5;
                e_step = dedx(proj, mat, c)*x_step*A_recip;
            }
            else if(x_step + x_traversed >= mat.thickness())
            {
                x_step = mat.thickness() - x_traversed;
                e_step = dedx(proj, mat, c)*x_step*A_recip;
                e_final -= e_step;
                proj.T = e_final;
                return (e_in - e_final)*proj.A;
            }
            else if(depth == s_depth_max)
            {
                return e_in*proj.A;
            }
            else if (e_final < 0.0)
            {
                return e_in * proj.A; //In case an integration step takes us below 0
            }
            else
            {
                e_step = dedx(proj, mat, c)*x_step*A_recip;
                e_final -= e_step;
                proj.T = e_final;
                x_traversed += x_step;
            }
        }
    }

    double reverse_integrate_energyloss(Projectile& proj, const Material& mat, const Config& c)
    {
        static double s_estep_max = 0.001;
        static int s_depth_max = 100;
        int depth = 0;
        double e_out = proj.T; //MeV/u
        double e_initial = e_out;
        double x_step = 0.25*mat.thickness(); //g/cm^2
        double x_traversed = 0.0;
        double e_step = dedx(proj, mat, c)*x_step;
		double A_recip = 1.0/proj.A;
        
        while(true)
        {
            if(e_step/e_initial > s_estep_max && depth < s_depth_max)
            {
                ++depth;
                x_step *= 0.5;
                e_step = dedx(proj, mat, c)*x_step*A_recip;
            }
            else if(x_step + x_traversed >= mat.thickness())
            {
                x_step = mat.thickness() -  x_traversed;
                e_step = dedx(proj, mat, c)*x_step;
                e_initial += e_step*A_recip;
                proj.T = e_initial;
                return (e_initial - e_out)*proj.A;
            }
            else if(depth == s_depth_max)
            {
                return e_out*proj.A;
            }
            else
            {
                e_step = dedx(proj, mat, c)*x_step;
                e_initial += e_step*A_recip;
                proj.T = e_initial;
                x_traversed += x_step;
            }
        }
    }
}
