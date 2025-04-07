#ifndef TOV_HPP_INCLUDED
#define TOV_HPP_INCLUDED

#include <vector>
#include <assert.h>
#include <optional>
#include <cmath>

double geometric_to_msol(double meters, double m_exponent);
double msol_to_geometric(double distance, double m_exponent);
double si_to_geometric(double quantity, double kg_exponent, double s_exponent);
double geometric_to_si(double quantity, double kg_exponent, double s_exponent);

///https://colab.research.google.com/drive/1yMD2j3Y6TcsykCI59YWiW9WAMW-SPf12#scrollTo=6vWjt7CWaVyV
///https://www.as.utexas.edu/astronomy/education/spring13/bromm/secure/TOC_Supplement.pdf
///https://arxiv.org/pdf/gr-qc/0403029

namespace tov
{
    template<typename T>
    double invert(T&& func, double y, double lower = 0, double upper = 1, bool should_search = true)
    {
        if(should_search)
        {
            //assumes upper > lower
            while(func(upper) < y)
                upper *= 2;
        }

        for(int i=0; i < 10000; i++)
        {
            double lower_mu = func(lower);
            double upper_mu = func(upper);

            double next = 0.5 * lower + 0.5 * upper;

            if(std::fabs(upper - lower) <= 1e-14)
                return next;

            //hit the limits of precision
            if(next == upper || next == lower)
                return next;

            double x = func(next);

            if(upper_mu >= lower_mu)
            {
                if(x >= y)
                    upper = next;
                ///x < y
                else
                    lower = next;
            }
            else
            {
                if(x >= y)
                    lower = next;
                else
                    upper = next;
            }
        }

        assert(false);
    };

    namespace eos
    {
        struct base
        {
            virtual double mu_to_p0(double mu) const{
                return invert(
                [&](double y){
                    return p0_to_mu(y);
                }, mu);
            };

            virtual double p0_to_mu(double p0) const = 0;

            virtual double mu_to_P(double mu) const{return p0_to_P(mu_to_p0(mu));};
            virtual double P_to_mu(double P) const{return p0_to_mu(P_to_p0(P));};

            virtual double P_to_p0(double P) const{
                return invert(
                [&](double y){
                    return p0_to_P(y);
                }, P);
            };

            virtual double p0_to_P(double p0) const = 0;
        };

        //units of c=g=msol=1
        struct polytrope : base
        {
            double Gamma = 0;
            double K = 0;

            polytrope(float _Gamma, float _K) : Gamma(_Gamma), K(_K) {}

            //virtual double mu_to_p0(double p0) const override;
            virtual double p0_to_mu(double p0) const override;

            //virtual double mu_to_P(double mu) const override;
            //virtual double P_to_mu(double P) const override;

            virtual double P_to_p0(double P) const override;
            virtual double p0_to_P(double p0) const override;
        };

        ///polytropic stuff
        ///https://www.as.utexas.edu/astronomy/education/spring13/bromm/secure/TOC_Supplement.pdf
        ///https://arxiv.org/pdf/0812.2163
        ///https://arxiv.org/pdf/2209.06052
    }

    struct integration_state
    {
        double m = 0;
        double p = 0;
    };

    integration_state make_integration_state(double central_rest_mass_density, double min_radius, const eos::base& param);
    integration_state make_integration_state_si(double central_rest_mass_density, double min_radius, const eos::base& param);

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

    std::optional<integration_solution> solve_tov(const integration_state& start, const eos::base& param, double min_radius, double min_pressure);
    ///returns a vector of central densities in units of c=g=msol, 1/length^2
    std::vector<double> search_for_rest_mass(double adm_mass, const eos::base& param);
}

namespace initial
{
    std::vector<double> calculate_isotropic_r(const tov::integration_solution& sol);
    std::vector<double> calculate_tov_phi(const tov::integration_solution& sol);
}

#endif // TOV_HPP_INCLUDED
