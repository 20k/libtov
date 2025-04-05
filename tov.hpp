#ifndef TOV_HPP_INCLUDED
#define TOV_HPP_INCLUDED

#include <vector>
#include <assert.h>
#include <functional>

double geometric_to_msol(double meters, double m_exponent);
double msol_to_geometric(double distance, double m_exponent);
double si_to_geometric(double quantity, double kg_exponent, double s_exponent);
double geometric_to_si(double quantity, double kg_exponent, double s_exponent);


///https://colab.research.google.com/drive/1yMD2j3Y6TcsykCI59YWiW9WAMW-SPf12#scrollTo=6vWjt7CWaVyV
///https://www.as.utexas.edu/astronomy/education/spring13/bromm/secure/TOC_Supplement.pdf
///https://arxiv.org/pdf/gr-qc/0403029

namespace tov
{
    namespace eos
    {
        struct numerical
        {
            ///linear mapping from 0 to max_mu
            ///mu is total specific energy density
            ///units of c=g=msol=1
            std::function<double(double)> p0_to_mu_func;
            std::function<double(double)> p0_to_P_func;

            double mu_to_p0(double p0) const;
            double p0_to_mu(double p0) const;

            double mu_to_P(double mu) const;
            double P_to_mu(double P) const;

            double P_to_p0(double P) const;
            double p0_to_P(double p0) const;
        };

        ///units of c=g=msol=1
        numerical from_polytropic(double Gamma, double K);
        ///https://www.as.utexas.edu/astronomy/education/spring13/bromm/secure/TOC_Supplement.pdf
        ///https://arxiv.org/pdf/0812.2163
        ///https://arxiv.org/pdf/2209.06052
        //numerical from_piecewise_polytropic();
    }

    struct integration_state
    {
        double m = 0;
        double p = 0;
    };

    integration_state make_integration_state(double central_rest_mass_density, double min_radius, const eos::numerical& param);
    integration_state make_integration_state_si(double central_rest_mass_density, double min_radius, const eos::numerical& param);

    struct integration_solution
    {
        double M_msol = 0;
        double R_msol = 0;

        std::vector<double> energy_density;
        std::vector<double> pressure;
        std::vector<double> cumulative_mass;
        std::vector<double> radius; //in schwarzschild coordinates, in units of c=G=mSol = 1

        int radius_to_index(double r) const;

        double M_geom() const;
        double R_geom() const;
    };

    integration_solution solve_tov(const integration_state& start, const eos::numerical& param, double min_radius, double min_pressure);
    ///returns a vector of central densities in units of c=g=msol, 1/length^2
    std::vector<double> search_for_rest_mass(double adm_mass, const eos::numerical& param);
}

namespace initial
{
    std::vector<double> calculate_isotropic_r(const tov::integration_solution& sol);
    std::vector<double> calculate_tov_phi(const tov::integration_solution& sol);
}

#endif // TOV_HPP_INCLUDED
