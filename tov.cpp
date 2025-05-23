#include "tov.hpp"
#include <cmath>
#include <numbers>
#include <stdio.h>
#include <iostream>

constexpr auto pi = std::numbers::pi;

static
double mix(double i1, double i2, double frac)
{
    return i1 * (1-frac) + i2 * frac;
}

///https://www.seas.upenn.edu/~amyers/NaturalUnits.pdf
//https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
double geometric_to_msol(double meters, double m_exponent)
{
    double m_to_kg = 1.3466 * std::pow(10., 27.);
    double msol_kg = 1.988416 * std::pow(10., 30.);
    double msol_meters = msol_kg / m_to_kg;

    return meters / pow(msol_meters, m_exponent);
}

double msol_to_geometric(double distance, double m_exponent)
{
    return geometric_to_msol(distance, -m_exponent);
}

double si_to_geometric(double quantity, double kg_exponent, double s_exponent)
{
    double G = 6.6743015 * pow(10., -11.);
    double C = 299792458;

    double factor = std::pow(G, -kg_exponent) * std::pow(C, 2 * kg_exponent - s_exponent);

    return quantity / factor;
}

double geometric_to_si(double quantity, double kg_exponent, double s_exponent)
{
    return si_to_geometric(quantity, -kg_exponent, -s_exponent);
}

double tov::eos::polytrope::p0_to_mu(double p0) const
{
    double P = p0_to_P(p0);

    return p0 + P / (Gamma - 1);
}

double tov::eos::polytrope::P_to_p0(double P) const
{
    return pow(P/K, 1/Gamma);
}

double tov::eos::polytrope::p0_to_P(double p0) const
{
    return K * pow(p0, Gamma);
}

template<typename T>
T interpolate_by_radius(const std::vector<double>& radius, const std::vector<T>& quantity, double r)
{
    assert(radius.size() >= 2);

    if(r <= radius.front())
        return quantity.front();

    if(r >= radius.back())
        return quantity.back();

    for(int i=0; i < (int)radius.size() - 1; i++)
    {
        double r1 = radius[i];
        double r2 = radius[i + 1];

        if(r > r1 && r <= r2)
        {
            double frac = (r - r1) / (r2 - r1);

            return mix(quantity[i], quantity[i + 1], frac);
        }
    }

    return quantity.back();
}

tov::integration_state tov::make_integration_state(double central_rest_mass_density, double rmin, const eos::base& param)
{
    double mu_c = param.p0_to_mu(central_rest_mass_density);

    double m = (4./3.) * pi * mu_c * std::pow(rmin, 3.);

    integration_state out;
    out.p = param.mu_to_P(mu_c);
    out.m = m;
    return out;
}

//p0 in si units
tov::integration_state tov::make_integration_state_si(double central_rest_mass_density, double rmin, const eos::base& param)
{
    //kg/m^3 -> m/m^3 -> 1/m^2
    double p0_geom = si_to_geometric(central_rest_mass_density, 1, 0);
    //m^-2 -> msol^-2
    double p0_msol = geometric_to_msol(p0_geom, -2);

    //std::cout << "density " << p0_msol << std::endl;

    return make_integration_state(p0_msol, rmin, param);
}

double tov::integration_solution::M_geom() const
{
    return msol_to_geometric(M_msol, 1);
}

double tov::integration_solution::R_geom() const
{
    return msol_to_geometric(R_msol, 1);
}

double tov::integration_solution::R_iso_msol() const
{
    return initial::calculate_isotropic_r(*this).back();
}

double tov::integration_solution::R_iso_geom() const
{
    return msol_to_geometric(R_iso_msol(), 1);
}

double tov::integration_solution::M0_msol() const
{
    double current = 0;
    double last_r = 0;

    for(int i=0; i < radius.size(); i++)
    {
        double r = radius[i];
        double dr = r - last_r;

        current += (4 * pi * r*r * rest_density[i] / sqrt(1 - 2 * cumulative_mass[i] / r)) * dr;

        last_r = r;
    }

    return current;
}

double tov::integration_solution::M0_geom() const
{
    return msol_to_geometric(M0_msol(), 1);
}

int tov::integration_solution::radius_to_index(double r) const
{
    assert(radius.size() > 0);

    if(r < radius.front())
        return 0;

    for(int i=1; i < radius.size(); i++)
    {
        if(r < radius[i])
            return i;
    }

    return radius.size() - 1;
}

struct integration_dr
{
    double dm = 0;
    double dp = 0;
};

std::optional<integration_dr> get_derivs(double r, const tov::integration_state& st, const tov::eos::base& param)
{
    std::optional<integration_dr> out;

    double mu = param.P_to_mu(st.p);

    double p = st.p;
    double m = st.m;

    //black hole + numerical error
    if(r <= 2 * m + 1e-12)
        return out;

    out.emplace();

    out->dm = 4 * pi * mu * std::pow(r, 2.);
    out->dp = -(mu + p) * (m + 4 * pi * r*r*r * p) / (r * (r - 2 * m));
    return out;
}

///units are c=g=msol=1
std::optional<tov::integration_solution> tov::solve_tov(const integration_state& start,  const tov::eos::base& param, double min_radius, double min_pressure)
{
    integration_state st = start;

    double current_r = min_radius;
    double dr = 1. / 1024.;

    integration_solution sol;

    double last_r = 0;
    double last_m = 0;

    while(1)
    {
        //save current integration state
        sol.rest_density.push_back(param.P_to_p0(st.p));
        sol.energy_density.push_back(param.P_to_mu(st.p));
        sol.pressure.push_back(st.p);
        sol.cumulative_mass.push_back(st.m);
        sol.radius.push_back(current_r);

        last_r = current_r;
        last_m = st.m;

        std::optional<integration_dr> data_opt = get_derivs(current_r, st, param);

        //oops, black hole!
        if(!data_opt.has_value())
            return std::nullopt;

        integration_dr& data = *data_opt;

        //euler integration
        st.m += data.dm * dr;
        st.p += data.dp * dr;
        current_r += dr;

        //something bad happened
        if(!std::isfinite(st.m) || !std::isfinite(st.p))
            return std::nullopt;

        //success!
        if(st.p <= min_pressure)
            break;
    }

    //sanity checks
    if(last_r <= min_radius * 100 || last_m < 0.0001f)
        return std::nullopt;

    sol.R_msol = last_r;
    sol.M_msol = last_m;

    return sol;
}

template<typename T>
inline
std::vector<double> search_for_mass(double mass, T&& get_mass)
{
    double rmin = 1e-6;
    double min_density = 0.;
    double max_density = 0.1;

    std::vector<double> densities;
    std::vector<std::optional<double>> masses;

    int to_check = 200;
    densities.resize(to_check);
    masses.resize(to_check);

    for(int i=0; i < to_check; i++)
    {
        double frac = (double)i / to_check;

        double test_density = mix(min_density, max_density, frac);

        std::optional<double> mass = get_mass(test_density);

        densities[i] = test_density;
        masses[i] = mass;
    }

    std::vector<double> out;

    for(int i=0; i < to_check - 1; i++)
    {
        if(!masses[i].has_value())
            continue;

        if(!masses[i+1].has_value())
            continue;

        double current = masses[i].value();
        double next = masses[i+1].value();

        double min_mass = std::min(current, next);
        double max_mass = std::max(current, next);

        if(mass >= min_mass && mass < max_mass)
        {
            auto get_cmass = [&](double in) {
                return get_mass(in).value();
            };

            double refined = tov::invert(get_cmass, mass, densities[i], densities[i + 1], false);

            out.push_back(refined);
        }
    }

    return out;
}

//personally i liked the voyage home better
std::vector<double> tov::search_for_adm_mass(double mass, const tov::eos::base& param)
{
    auto get_mass = [&](double density_in) -> std::optional<double> {
        integration_state st = make_integration_state(density_in, 1e-6, param);
        auto next_sol_opt = solve_tov(st, param, 1e-6, 0.);

        std::optional<double> ret;

        if(next_sol_opt)
            ret = next_sol_opt->M_msol;

        return ret;
    };

    return search_for_mass(mass, get_mass);
}

std::vector<double> tov::search_for_rest_mass(double mass, const tov::eos::base& param)
{
    auto get_mass = [&](double density_in) -> std::optional<double> {
        integration_state st = make_integration_state(density_in, 1e-6, param);
        auto next_sol_opt = solve_tov(st, param, 1e-6, 0.);

        std::optional<double> ret;

        if(next_sol_opt)
            ret = next_sol_opt->M0_msol();

        return ret;
    };

    return search_for_mass(mass, get_mass);
}

std::vector<double> initial::calculate_isotropic_r(const tov::integration_solution& sol)
{
    std::vector<double> r_hat;
    double last_r = 0;
    double log_rhat_r = 0;

    for(int i=0; i < (int)sol.radius.size(); i++)
    {
        double r = sol.radius[i];
        double m = sol.cumulative_mass[i];

        double dr = r - last_r;

        double rhs = (1 - sqrt(1 - 2 * m/r)) / (r * sqrt(1 - 2 * m/r));

        log_rhat_r += dr * rhs;

        double lr_hat = exp(log_rhat_r);

        r_hat.push_back(lr_hat);

        last_r = r;
    }

    {
        double final_r = r_hat.back();

        double R = sol.radius.back();
        double M = sol.cumulative_mass.back();

        double scale = (1/(2*R)) * (sqrt(R*R - 2 * M * R) + R - M) / final_r;

        for(int i=0; i < (int)sol.radius.size(); i++)
        {
            r_hat[i] *= sol.radius[i] * scale;
        }
    }

    return r_hat;
}

std::vector<double> initial::calculate_tov_phi(const tov::integration_solution& sol)
{
    auto isotropic_r = calculate_isotropic_r(sol);
    int samples = sol.radius.size();

    std::vector<double> phi;
    phi.resize(samples);

    for(int i=0; i < (int)sol.radius.size(); i++)
    {
        phi[i] = pow(sol.radius[i] / isotropic_r[i], 1./2.);
    }

    #if 0
    auto isotropic_to_schwarzschild = [&](auto isotropic_in)
    {
        return interpolate_by_radius(isotropic_r, sol.radius, isotropic_in);
    };

    double check_mass = 0;
    double last_r_bar = 0;

    for(int i=0; i < (int)samples; i++)
    {
        double r_bar = isotropic_r[i];

        double r = isotropic_to_schwarzschild(r_bar);
        double e = tov::interpolate_by_radius(sol.radius, sol.energy_density, r);

        check_mass += 4 * pi * r_bar*r_bar * pow(phi[i], 5.) * e * (r_bar - last_r_bar);

        last_r_bar = r_bar;
    }

    std::cout << "Check Mass " << check_mass << std::endl;
    #endif

    return phi;
}
