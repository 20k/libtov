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
    tov::eos::polytrope param(2, 123.641);

    double rmin = 1e-6;

    double min_density = 0;
    double max_density = 0.1;
    int to_check = 20000;

    std::string str = "Mass, Radius\n";

    for(int i=0; i < to_check; i++)
    {
        double frac = (double)i / to_check;

        double test_density = mix(min_density, max_density, frac);

        tov::integration_state next_st = tov::make_integration_state(test_density, rmin, param);
        auto next_sol_opt = tov::solve_tov(next_st, param, rmin, 0.);

        if(!next_sol_opt)
            continue;

        tov::integration_solution& next_sol = *next_sol_opt;

        str += std::to_string(next_sol.M_msol) + ", " + std::to_string(next_sol.R_geom()/1000.) + "\n";

        //std::cout << next_sol.M_msol << ", " << next_sol.R_geom() / 1000. << std::endl;
    }

    std::ofstream out("data.csv");
    out << str;
}

void test_1()
{
    tov::eos::polytrope param(2, 123.641);
    double paper_p0 = 6.235 * pow(10., 17.);

    double rmin = 1e-6;

    tov::integration_state st = tov::make_integration_state_si(paper_p0, rmin, param);

    tov::integration_solution sol = tov::solve_tov(st, param, rmin, 0).value();

    std::cout << "Solved for " << sol.R_geom() / 1000. << "km " << sol.M_msol << " msols " << std::endl;
}

void test_2()
{
    ///should find 0.00100957
    tov::eos::polytrope param(2, 123.641);

    auto results = tov::search_for_adm_mass(1.543, param);

    for(auto& i : results)
    {
        std::cout << "Density " << i << std::endl;
    }
}

///https://arxiv.org/pdf/gr-qc/0110047
void test_3()
{
    tov::eos::polytrope param(2, 100);
    double paper_p0 = 8 * pow(10., -3.);

    double rmin = 1e-6;

    tov::integration_state st = tov::make_integration_state(paper_p0, rmin, param);

    tov::integration_solution sol = tov::solve_tov(st, param, rmin, 0).value();

    std::cout << "Solved for " << sol.R_geom() / 1000. << "km " << sol.M_msol << " msols m0: " << sol.M0_msol() << " msols " << std::endl;
}

int main()
{
    //mass_radius_curve();

    test_1();
    test_3();

    test_2();

    //kg/m^3

    //std::vector<double> tov_phi = initial::calculate_tov_phi(sol);

    return 0;
}
