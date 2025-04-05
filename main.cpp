#include <fstream>
#include <string>
#include <iostream>
#include "tov.hpp"
#include <cmath>
#include <numbers>

constexpr auto pi = std::numbers::pi;

static
double mix(double i1, double i2, double frac)
{
    return i1 * (1-frac) + i2 * frac;
}

void mass_radius_curve()
{
    /*tov::parameters param;
    param.K = 123.641;
    param.Gamma = 2;*/

    tov::eos::numerical param = tov::eos::from_polytropic(2, 123.641);

    double adm_mass = 1.5;

    double r_approx = adm_mass / 0.06;

    double start_E = adm_mass / ((4./3.) * pi * r_approx*r_approx*r_approx);
    double start_P = param.mu_to_P(start_E);
    double start_density = param.P_to_p0(start_P);

    double rmin = 1e-6;

    double min_density = start_density / 100;
    double max_density = start_density * 2000;
    int to_check = 20000;

    std::string str = "Mass, Radius\n";

    for(int i=0; i < to_check; i++)
    {
        double frac = (double)i / to_check;

        double test_density = mix(min_density, max_density, frac);

        tov::integration_state next_st = tov::make_integration_state(test_density, rmin, param);
        tov::integration_solution next_sol = tov::solve_tov(next_st, param, rmin, 0.);

        str += std::to_string(next_sol.M_msol) + ", " + std::to_string(next_sol.R_geom()/1000.) + "\n";

        //std::cout << next_sol.M_msol << ", " << next_sol.R_geom() / 1000. << std::endl;
    }

    std::ofstream out("data.csv");
    out << str;
}

int main()
{
    tov::eos::numerical param = tov::eos::from_polytropic(2, 123.641);

    /*return;
    auto results = tov::search_for_adm_mass(1.543, param);

    for(auto& i : results)
    {
        std::cout << "Density " << i << std::endl;
    }

    assert(false);*/

    //kg/m^3

    double paper_p0 = 6.235 * pow(10., 17.);

    //this is in c=g=msol, so you'd need to use make_integration_state()
    //double p0 = 1.28e-3;


    double rmin = 1e-6;

    //integration_state st = make_integration_state(p0, rmin, param);
    tov::integration_state st = tov::make_integration_state_si(paper_p0, rmin, param);

    tov::integration_solution sol = tov::solve_tov(st, param, rmin, 0);

    std::cout << "Solved for " << sol.R_geom() / 1000. << "km " << sol.M_msol << " msols " << std::endl;

    std::vector<double> tov_phi = initial::calculate_tov_phi(sol);

    return 0;
}
