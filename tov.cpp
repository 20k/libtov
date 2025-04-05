#include "tov.hpp"
#include <cmath>
#include <numbers>

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

double table_lookup(double val, const std::vector<double>& val_is_in, const std::vector<double>& companion)
{
    assert(val >= 0);
    assert(val_is_in.size() > 0);
    assert(companion.size() > 0);
    assert(val_is_in.size() == companion.size());

    assert(val >= val_is_in.front());
    assert(val <= val_is_in.back());

    if(val == val_is_in.front())
        return companion.front();

    if(val == val_is_in.back())
        return companion.back();

    for(int i=0; i < (int)companion.size() - 1; i++)
    {
        double v_1 = val_is_in[i];
        double v_2 = val_is_in[i + 1];

        if(val > v_1 && val <= v_2)
        {
            double frac = (val - v_1) / (v_2 - v_1);

            double o_1 = companion[i];
            double o_2 = companion[i + 1];

            return mix(o_1, o_2, frac);
        }
    }

    assert(false);
}

double invert(std::function<double(double)> func, double y)
{
    double lower = 0;
    double upper = 1;

    while(func(upper) < y)
        upper *= 2;

    for(int i=0; i < 1024; i++)
    {
        double lower_mu = func(lower);
        double upper_mu = func(upper);

        double next = (lower + upper)/2.;

        double x = func(next);

        if(x >= y)
        {
            upper = next;
        }
        ///x < y
        else
        {
            lower = next;
        }
    }

    return (lower + upper)/2;
};


double tov::eos::numerical::mu_to_p0(double mu) const
{
    return invert(p0_to_mu_func, mu);
}

double tov::eos::numerical::p0_to_mu(double p0) const
{
    return p0_to_mu_func(p0);
}

double tov::eos::numerical::mu_to_P(double mu) const
{
    double p0 = invert(p0_to_mu_func, mu);

    return p0_to_P_func(p0);
}

double tov::eos::numerical::P_to_mu(double P) const
{
    double p0 = invert(p0_to_P_func, P);
    return p0_to_mu_func(p0);
}

double tov::eos::numerical::P_to_p0(double P) const
{
    return invert(p0_to_P_func, P);
}

double tov::eos::numerical::p0_to_P(double p0) const
{
    return p0_to_P_func(p0);
}

tov::eos::numerical tov::eos::from_polytropic(double Gamma, double K)
{
    tov::eos::numerical out;

    auto rest_mass_density_to_pressure = [=](double rest_mass_density)
    {
        return K * pow(rest_mass_density, Gamma);
    };

    auto rest_mass_density_to_energy_density = [=](double rest_mass_density)
    {
        double p = rest_mass_density_to_pressure(rest_mass_density);

        double p0 = rest_mass_density;

        return p0 + p/(Gamma-1);
    };

    out.p0_to_mu_func = rest_mass_density_to_energy_density;
    out.p0_to_P_func = rest_mass_density_to_pressure;

    return out;
}

template<typename T>
inline
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

#if 0
double tov::parameters::rest_mass_density_to_pressure(double rest_mass_density) const
{
    return K * pow(rest_mass_density, Gamma);
}

double tov::parameters::rest_mass_density_to_energy_density(double rest_mass_density) const
{
    double p = rest_mass_density_to_pressure(rest_mass_density);

    double p0 = rest_mass_density;

    return p0 + p/(Gamma-1);
}

///inverse equation of state
///p -> p0
double tov::parameters::pressure_to_rest_mass_density(double p) const
{
    return std::pow(p/K, 1/Gamma);
}

///e = p0 + P/(Gamma-1)
double tov::parameters::pressure_to_energy_density(double p) const
{
    return pressure_to_rest_mass_density(p) + p / (Gamma - 1);
}

double tov::parameters::energy_density_to_pressure(double mu) const
{
    auto func = [&](double arg)
    {
        return rest_mass_density_to_energy_density(arg);
    };

    ///lets solve this numerically
    ///mu = p0 + P/(Gamma-1)
    double lower = 0;
    double upper = 1;

    while(func(upper) < mu)
        upper *= 2;

    for(int i=0; i < 1024; i++)
    {
        double lower_mu = func(lower);
        double upper_mu = func(upper);

        double next = (lower + upper)/2.;

        double next_mu = func(next);

        if(next_mu >= mu)
        {
            upper = next;
        }
        ///next_mu < mu
        else
        {
            lower = next;
        }
    }

    return (lower + upper)/2;
}
#endif

tov::integration_state tov::make_integration_state(double central_rest_mass_density, double rmin, const eos::numerical& param)
{
    double mu_c = param.p0_to_mu(central_rest_mass_density);

    double m = (4./3.) * pi * mu_c * std::pow(rmin, 3.);

    integration_state out;
    out.p = param.mu_to_P(mu_c);
    out.m = m;
    return out;
}

//p0 in si units
tov::integration_state tov::make_integration_state_si(double central_rest_mass_density, double rmin, const eos::numerical& param)
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

integration_dr get_derivs(double r, const tov::integration_state& st, const tov::eos::numerical& param)
{
    double e = param.P_to_mu(st.p);

    double p = st.p;
    double m = st.m;

    integration_dr out;

    out.dm = 4 * pi * e * std::pow(r, 2.);
    out.dp = -(e + p) * (m + 4 * pi * r*r*r * p) / (r * (r - 2 * m));
    return out;
}

///units are c=g=msol
tov::integration_solution tov::solve_tov(const integration_state& start,  const tov::eos::numerical& param, double min_radius, double min_pressure)
{
    integration_state st = start;

    double current_r = min_radius;
    double dr = 1. / 1024.;

    integration_solution sol;

    double last_r = 0;
    double last_m = 0;

    while(1)
    {
        sol.energy_density.push_back(param.P_to_mu(st.p));
        sol.pressure.push_back(st.p);
        sol.cumulative_mass.push_back(st.m);

        double r = current_r;

        sol.radius.push_back(r);

        last_r = r;
        last_m = st.m;

        integration_dr data = get_derivs(r, st, param);

        st.m += data.dm * dr;
        st.p += data.dp * dr;
        current_r += dr;

        if(!std::isfinite(st.m) || !std::isfinite(st.p))
            break;

        if(st.p <= min_pressure)
            break;
    }

    sol.R_msol = last_r;
    sol.M_msol = last_m;

    return sol;
}

//personally i liked the voyage home better
std::vector<double> tov::search_for_rest_mass(double adm_mass, const tov::eos::numerical& param)
{
    double r_approx = adm_mass / 0.06;

    double start_E = adm_mass / ((4./3.) * pi * r_approx*r_approx*r_approx);
    double start_P = param.mu_to_P(start_E);
    double start_density = param.mu_to_p0(param.P_to_mu(start_P));

    double rmin = 1e-6;

    std::vector<double> densities;
    std::vector<double> masses;

    int to_check = 2000;
    densities.resize(to_check);
    masses.resize(to_check);

    double min_density = start_density / 100;
    double max_density = start_density * 200;

    for(int i=0; i < to_check; i++)
    {
        double frac = (double)i / to_check;

        double test_density = mix(min_density, max_density, frac);

        integration_state next_st = make_integration_state(test_density, rmin, param);
        integration_solution next_sol = solve_tov(next_st, param, rmin, 0.);

        densities[i] = test_density;
        masses[i] = next_sol.M_msol;
    }

    std::vector<double> out;

    for(int i=0; i < to_check - 1; i++)
    {
        double current = masses[i];
        double next = masses[i+1];

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
