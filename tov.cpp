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
    double e = param.P_to_mu(st.p);

    double p = st.p;
    double m = st.m;

    //black hole + numerical error
    if(r <= 2 * m + 1e-12)
        return std::nullopt;

    integration_dr out;

    out.dm = 4 * pi * e * std::pow(r, 2.);
    out.dp = -(e + p) * (m + 4 * pi * r*r*r * p) / (r * (r - 2 * m));
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

    //int i = 0;

    while(1)
    {
        sol.energy_density.push_back(param.P_to_mu(st.p));
        sol.pressure.push_back(st.p);
        sol.cumulative_mass.push_back(st.m);

        double r = current_r;

        sol.radius.push_back(r);

        /*if((i % 1000) == 0)
        //if(i > 10000)
            printf("R %.24f m %.24f\n", r, st.m);

        i++;*/

        last_r = r;
        last_m = st.m;

        std::optional<integration_dr> data_opt = get_derivs(r, st, param);

        if(!data_opt.has_value())
            return std::nullopt;

        integration_dr& data = *data_opt;

        st.m += data.dm * dr;
        st.p += data.dp * dr;
        current_r += dr;

        if(!std::isfinite(st.m) || !std::isfinite(st.p))
            return std::nullopt;

        if(st.p <= min_pressure)
        {
            //printf("Pressure broke\n");
            break;
        }
    }

    //printf("Fr %.24f %.24f min %.24f\n", last_r, last_m, min_radius);

    if(last_r <= min_radius * 100 || last_m < 0.0001f)
    {
        //printf("Nope\n");
        return std::nullopt;
    }

    sol.R_msol = last_r;
    sol.M_msol = last_m;

    return sol;
}

//personally i liked the voyage home better
std::vector<double> tov::search_for_rest_mass(double adm_mass, const tov::eos::base& param)
{
    //todo: this function is a mess
    //lets imagine that all the mass were crammed into the schwarzschild radius
    //we have a uniform density. This is the minimum radius of the neutron star
    //so, as pc goes up, the mass : radius ratio will start to approach a black hole
    //we can use this ratio as our upper bound
    //schwarzschild radius = 2M
    //the issue is, we can't really a priori know if something's a black hole
    //need to cooperatively solve with the tov solver, as the bottom term (r * (r - 2 * m))  denotes it being a black hole
    /*double r_approx = adm_mass / 0.06;

    double start_E = adm_mass / ((4./3.) * pi * r_approx*r_approx*r_approx);
    double start_P = param.mu_to_P(start_E);
    double start_density = param.mu_to_p0(param.P_to_mu(start_P));*/

    /*double rmin = 1e-6;

    double min_density = 0;
    double max_density = 0;

    {
        double r_approx = adm_mass / 0.06;

        double start_E = adm_mass / ((4./3.) * pi * r_approx*r_approx*r_approx);
        double start_P = param.mu_to_P(start_E);
        double start_density = param.mu_to_p0(param.P_to_mu(start_P));

        double d1 = start_density;
        double d2 = start_density;

        //std::cout << "d1 " << d1 << std::endl;

        tov::integration_state st1 = make_integration_state(d1, rmin, param);

        tov::integration_solution sol1 = tov::solve_tov(st1, param, rmin, 0.f).value();

        tov::integration_state st2 = make_integration_state(d2, rmin, param);
        std::optional<tov::integration_solution> sol2 = tov::solve_tov(st2, param, rmin, 0.f);

        while(sol2)
        {
            d2 *= 2;

            printf("Testing density %.24f\n", d2);

            st2 = make_integration_state(d2, rmin, param);
            sol2 = tov::solve_tov(st2, param, rmin, 0.f);
        }

        printf("d Min %.24f max %.24f\n", d1, d2);

        //assert(false);

        ///in theory, d1 is valid, d2 is schwarzschild. We want to find the exact bound which is valid, and use that as max density

        for(int i=0; i < 128; i++)
        {
            double next = (d1 + d2) / 2.;

            tov::integration_state st = make_integration_state(next, rmin, param);
            auto tov_next = tov::solve_tov(st, param, rmin, 0.f);

            if(tov_next)
            {
                d1 = next;
            }
            else
            {
                d2 = next;
            }
        }

        min_density = 0;
        max_density = d2;
    }*/

    double rmin = 1e-6;
    double min_density = 0.f;
    double max_density = 1.f;

    printf("Min %.24f max %.24f\n", min_density, max_density);

    //assert(false);

    std::vector<double> densities;
    std::vector<std::optional<double>> masses;

    int to_check = 2000;
    densities.resize(to_check);
    masses.resize(to_check);

    for(int i=0; i < to_check; i++)
    {
        printf("Checking %i\n", i);

        double frac = (double)i / to_check;

        double test_density = mix(min_density, max_density, frac);

        integration_state next_st = make_integration_state(test_density, rmin, param);
        std::optional<integration_solution> next_sol_opt = solve_tov(next_st, param, rmin, 0.);

        std::optional<double> mass;

        if(next_sol_opt)
            mass = next_sol_opt->M_msol;

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

        if(adm_mass >= min_mass && adm_mass < max_mass)
        {
            double frac = (adm_mass - min_mass) / (max_mass - min_mass);

            out.push_back(mix(densities[i], densities[i+1], frac));
        }
    }

    return out;
}

std::vector<double> initial::calculate_isotropic_r(const tov::integration_solution& sol)
{
    std::vector<double> dlog_dr;
    dlog_dr.reserve(sol.cumulative_mass.size());

    for(int i=0; i < (int)sol.cumulative_mass.size(); i++)
    {
        double r = sol.radius[i];
        double m = sol.cumulative_mass[i];

        double rhs = (1 - sqrtf(1 - 2 * m/r)) / (r * sqrt(1 - 2 * m/r));

        //double rhs = (pow(r, 0.5) - pow(r - 2 * m, 0.5)) / (r * pow(r - 2 * m, 0.5));
        dlog_dr.push_back(rhs);
    }

    std::vector<double> r_hat;
    double last_r = 0;
    double log_rhat_r = 0;

    for(int i=0; i < (int)sol.radius.size(); i++)
    {
        double r = sol.radius[i];

        double dr = r - last_r;

        log_rhat_r += dr * dlog_dr[i];

        //std::cout << "step size " << dr * dlog_dr[i] << std::endl;

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

    auto isotropic_to_schwarzschild = [&](auto isotropic_in)
    {
        return interpolate_by_radius(isotropic_r, sol.radius, isotropic_in);
    };

    int samples = sol.radius.size();

    std::vector<double> phi;
    phi.resize(samples);

    for(int i=0; i < (int)sol.radius.size(); i++)
    {
        phi[i] = pow(sol.radius[i] / isotropic_r[i], 1./2.);
    }

    #if 0
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
