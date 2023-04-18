/*
 * Testing for gwm_integrators
 *
 * Dec 2022, Gordon McCann
 */
#include "catima/gwm_integrators.h"
#include "catima/nucdata.h"
#include <iostream>

int main()
{
    std::cout<<"-------Testing GWM Energy Loss Integration-------"<<std::endl;
    catima::Projectile p1(catima::element_atomic_weight(1), 1.0, 0, 3.0);
    catima::Material mat1(catima::get_material(6));
    mat1.density(2.23).thickness(500.0*1e-6);

    double result = catima::integrate_energyloss(p1, mat1);
    std::cout<<"Energy loss (MeV): "<<result<<" Final energy: "<<p1.T<<std::endl;
    result = catima::reverse_integrate_energyloss(p1, mat1);
    std::cout<<"Reverse Energy loss (MeV): "<<result<<" Initial energy: "<<p1.T<<std::endl;

}
