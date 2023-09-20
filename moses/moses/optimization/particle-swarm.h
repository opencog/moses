/*
 * moses/moses/optimization/particle-swarm.h
 *
 * Copyright (C) 2002-2008 Novamente LLC
 * Copyright (C) 2012 Poulin Holdings
 * All Rights Reserved
 *
 * Written by Arley Ristar
 *            Nil Geisweiller
 *            Linas Vepstas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License v3 as
 * published by the Free Software Foundation and including the exceptions
 * at http://opencog.org/wiki/Licenses
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program; if not, write to:
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef _MOSES_PARTICLE_SWARM_H
#define _MOSES_PARTICLE_SWARM_H

#include <opencog/util/oc_assert.h>

#include "../representation/instance_set.h"
#include "optimization.h"

namespace opencog { namespace moses {

typedef std::vector<double> velocity;

// Half of the things implemented here (about PSO) is based on or citatios of this article:
// Nouaouria, Nabila, and Mounir Boukadoum. "Improved global-best particle swarm
// optimization algorithm with mixed-attribute data classification capability."
// Applied Soft Computing 21 (2014): 554-567.

// Particle Swarm parameters
// XXX PSO parameters hardcoded just for now.
struct ps_parameters
{
    // There isn't need to set all this parameters, for most
    // problems it works with the default values.
    ps_parameters()
        : max_parts(50),
          bit_c1(0.7),
          disc_c1(2.05),
          cont_c1(0.7),
          bit_c2(1.43),
          disc_c2(2.05),
          cont_c2(1.43),
          inertia(0.7), // static inertia
          bit_min_value(0),
          bit_max_value(1),
          disc_min_value(0),
          disc_max_value(1),
          cont_min_value(std::numeric_limits<int>::min()),
          cont_max_value(std::numeric_limits<int>::max()),
          bit_min_vel(-6), // 1.8%
          bit_max_vel(6), // 98.2%
          disc_min_vel(-disc_max_value/2),
          disc_max_vel(disc_max_value/2),
          cont_min_vel(std::numeric_limits<int>::min() / 2),
          cont_max_vel(std::numeric_limits<int>::max() / 2),
          range_bit_vel(bit_max_vel - bit_min_vel),
          range_cont_vel(cont_max_vel - cont_min_vel)
    {
        double fi = disc_c1 + disc_c2; // where c1 + c2 > 4
        constriction_disc = 2 / (2 - fi - std::sqrt((fi * fi) - (4 * fi)));
    }

    // Maximum number of particles per deme
    unsigned max_parts;

    // For continuous problems, 0.7 and 1.43 are good values.
    // XXX I have to find the best values for bit and disc.
    // Information from: M. Clerc, L’optimisation par essaim particulaire: versions
    // paramétriques et adaptatives, Hermes Science Publications, Lavoisier, Paris, 2005.
    double bit_c1, disc_c1, cont_c1; // c1 = Individual learning rate.
    double bit_c2, disc_c2, cont_c2; // c2 = Social parameter.

    // See Inertia Weight Strategies in Particle Swarm Optimization
    // For better ways than static.
    // Inertia is just used in continuous for now.
    double inertia; // Static Inertia for now.

    // Values and Velocity:
    // For bit:
    // Use DPSO (original) where the bit have a minimal value of 0 and maximal of 1.
    // The concept of velocity doesn't exist as in continuous. Now it represents
    // changes of probabilities that a bit will be in one state or the other.
    // From here:
    // J. Kennedy, R. Eberhart, A discrete binary version of the particle swarm
    // algorithm, IEEE Conf. Syst. Man Cybern. 5 (1997) 4104–4108.
    //
    // For discrete:
    // Use rounding off for now XXX it isn't so effective, but for now is easier
    // to implement. I'll use MVPSO, JPSO or IMVPSO later.
    // In the transformation i'll use min value of 0 and 1 for max in the continuous
    // space mapped from 0 to .multy() in discrete space.
    // For velocity, i'll use the commom rule [-range/2, range/2].
    // In this case: [-0.5, 0.5]
    // For rounding off, i'm using this modification:
    // H.A. Hassan, I.M. Yassin, A.K. Halim, A. Zabidi, Z.A. Majid, H.Z. Abidin, Logical
    // effort using a novel discrete particle swarm optimization algorithm, in: 5th
    // International Colloquium on Signal Processing & Its Applications (CSPA), 2009,
    // 978-1-4244-4152-5/09.
    //
    // For continuous:
    // Classical PSO from here:
    // J. Kennedy, R. Eberhart, Particle swarm optimization, in: Proceedings of the 4th
    // IEEE International Conference on Neural Networks, Perth, Australia, 1995, pp.
    // 1942–1948.
    // With two more mechanism: confinement and wind dispersion from:
    // Confinement:
    // M. Clerc, L’optimisation par essaim particulaire: versions paramétriques et
    // adaptatives, Hermes Science Publications, Lavoisier, Paris, 2005.
    // Wind Dispersion:
    // K. Chandramouli, E. Izquierdo, Image classification using chaotic particle swarm
    // optimization, in: Proceedings of the International Conference on Image Pro-
    // cessing (ICIP ‘06), 2006.
    // The two combined:
    // N. Nouaouria, M. Boukadoum, Particle swarm classification for high dimen-
    // sional data sets, in: Proceedings of 22th International IEEE Conference on Tools
    // with Artificial Intelligence IEEE-ICTAI, vol. 1, Arras, France, 27–29 October 2010,
    // 2010, pp. 87–93.
    // XXX If i have time, i'll put some variation here to get a better global search.
    double bit_min_value, bit_max_value, // [0,1] <- XXX these two aren't used yet.
           disc_min_value, disc_max_value, // [0,1] in rounding off
           cont_min_value, cont_max_value; // [min contin_t, max contin_t]
    double bit_min_vel, bit_max_vel, // [0,1]
           disc_min_vel, disc_max_vel, // [-0.5,0.5] in rounding off
           cont_min_vel, cont_max_vel; // [min contin_t, max contin_t]
    double range_bit_vel, range_cont_vel;

    // From:
    // H.A. Hassan, I.M. Yassin, A.K. Halim, A. Zabidi, Z.A. Majid, H.Z. Abidin, Logical
    // effort using a novel discrete particle swarm optimization algorithm, in: 5th
    // International Colloquium on Signal Processing & Its Applications (CSPA), 2009,
    // 978-1-4244-4152-5/09
    double constriction_disc;
};

////////////////////
// Particle Swarm //
////////////////////

// TODO: pso description
struct particle_swarm : optimizer_base
{
    particle_swarm(const optim_parameters& op = optim_parameters(),
                    const ps_parameters& ps = ps_parameters())
        : optimizer_base(op), _total_RAM_bytes(getTotalRAM()), ps_params(ps) {}

protected:
    // Variables:
    const uint64_t _total_RAM_bytes;
    size_t _instance_bytes;
    const ps_parameters ps_params;

    // Structs
    // Discrete values have to be outside the deme.
    struct discrete_particles {
        std::vector<std::vector<double>> temp, best_personal;
        discrete_particles(unsigned part_size, unsigned disc_size) :
            temp(part_size, std::vector<double>(disc_size)),
            best_personal(part_size, std::vector<double>(disc_size)) {}
    };

    // Functions (Better explanation in declaration):
    // log legend for graph stats
    void log_stats_legend();

    unsigned calc_swarm_size (const field_set& fs);

//    void brute_force(deme_t deme, const instance& init_inst,
//        const field_set& fs, const iscorer_base& iscorer);

    void initialize_particles (const unsigned& swarm_size,
        deme_t& best_parts, std::vector<velocity>& velocities,
        discrete_particles& disc_parts, const field_set& fields);

    void initialize_random_particle (instance& new_inst, velocity& vel,
        std::vector<double>& dist_values, const field_set& fs);

////// Velocity Functions //////
    //// Check bounds functions:
    // There's no real bounds check for bit velocity, it's the probability.
    void check_bit_vel(double &vel) { // Check bounds of bit velocity
        vel = clamp(vel, ps_params.bit_min_vel, ps_params.bit_max_vel); }
    void check_disc_vel(double &vel) { // Check bounds of a discrete velocity
        vel = clamp(vel, ps_params.disc_min_vel, ps_params.disc_max_vel); }
    void check_cont_vel(double &vel) { // Check bounds of a continuous velocity
        vel = clamp(vel, ps_params.cont_min_vel, ps_params.cont_max_vel); }

    //// Generate initial random velocity
    double gen_bit_vel() { // [0,1]
        return (randGen().randdouble() *
                ps_params.range_bit_vel) - ps_params.bit_max_vel;
    }
    double gen_disc_vel() { // [-0.5, 0.5] when mapped
        return (randGen().randdouble() - ps_params.disc_max_vel);
    }
    double gen_cont_vel() {
        return (randGen().randdouble() *
                ps_params.range_cont_vel ) - ps_params.cont_max_vel;
    }

////// Particle values functions //////
    //// Generate initial random instance knob value
    bool gen_bit_value() { // {0,1}
        return randGen().randbool(); }
    double gen_disc_value() { // [0,1]
        return randGen().randdouble(); }
    double gen_cont_value() { //
        return (randGen().randdouble() *
                ps_params.range_cont_vel) + ps_params.cont_min_vel;
    }

    //// Confinament functions:
    // There isn't confinement for bit values
    // The update rule already kind of do it.
    void confinement_disc(double& value) {
        value = clamp(value, ps_params.disc_min_value, ps_params.disc_max_value);
    }

    void confinement_cont(double& value) {
        value = clamp(value, ps_params.cont_min_value, ps_params.cont_max_value);
    }

////// Update specific functions //////
    //// All
    void update_particles(deme_t& temp_parts, const deme_t& best_parts, const int& best_index,
        std::vector<velocity>& velocities, discrete_particles& disc_parts, const field_set& fields);

    //// Bit
    //
    void update_bit_vel(double& vel, int&& temp,
            int&& personal, int&& global) { // Bit type is bool
        // Bool convertion to int: false to 0, true to 1.
        // Bit vel hasn't inertia.
        vel += (randGen().randdouble() * (personal - temp)) +
            (randGen().randdouble() * (global - temp));
        check_bit_vel(vel);
    }

    // XXX Explanation
    bool new_bit_value(const double& vel){
        return (randGen().randdouble() < // Sigmoid
                (1 / (1 + std::exp(-vel)))); // XXX if slow try f(x) = x / (1 + abs(x)) or tanh(x)
    }

    void update_bit_particle(instance& temp, const instance& personal,
        const instance& global, velocity::iterator vel, const field_set& fs);

    //// Discrete
    // Update discrete velocity
    void update_disc_vel(double& vel, const double& temp,
            const double& personal, const double& global) { // Disc type before conversion is double
        // Disc vel hasn't inertia, but has Constriction Factor
        vel += (ps_params.disc_c1 * randGen().randdouble() * (personal - temp)) +
            (ps_params.disc_c2 * randGen().randdouble() * (global - temp));
        vel = vel * ps_params.constriction_disc; // Constriction part
        check_disc_vel(vel);
    }

    disc_t cont2disc(double& cvalue, const unsigned max_dvalue){
        // The original formula is dvalue = std::round(
        // (((max_dvalue - min_dvalue) * (cvalue - min_cvalue))/(max_cvalue - min_cvalue)) + min_dvalue);
        // But [min_cvalue, max_cvalue] == [0,1] and
        // min_dvalue == 0 and max_dvalue == it.multy(), so...
        return (disc_t) std::round(cvalue * (max_dvalue - 1)); // Return dvalue
    }

    // XXX Explanation
    disc_t new_disc_value(double& cvalue,
            const double& vel, const unsigned max_dvalue){
        cvalue += vel;
        confinement_disc(cvalue);
        return cont2disc(cvalue, max_dvalue);
    }

    // Update discrete part of the particle
    void update_disc_particle(instance& dtemp, std::vector<double>& temp,
        const std::vector<double>& personal, const std::vector<double>& global,
        velocity::iterator vel, const field_set& fs);

    //// Continuous
    // Update continuous velocity
    void update_cont_vel(double& vel, const double& temp,
            const double& personal, const double& global) { // Cont type is double
        vel += ps_params.inertia * ((ps_params.cont_c1 * randGen().randdouble() * (personal - temp)) +
                (ps_params.cont_c2 * randGen().randdouble() * (global - temp)));
        check_cont_vel(vel);
    }

    // XXX Explanation
    contin_t new_cont_value(const contin_t& value, const double& vel){
        // Wind dispersion enter here.
        contin_t res = value + vel;
        confinement_cont(res);
        return res;
    }

    // Update contin part of the particle
    void update_cont_particle(instance& temp, const instance& personal,
            const instance& global, velocity::iterator vel, const field_set& fs);

    // TODO: Wind dispersion, but test without first
    // Make it later is easy.

public:
    /**
     * XXX Perform search of the local neighborhood of an instance.  The
     * search is exhaustive if the local neighborhood is small; else
     * the local neighborhood is randomly sampled.
     *
     * @param deme      Where to store the candidates searched. The deme
     *                  is assumed to be empty.  If it is not empty, it
     *                  will be overwritten.
     * @prama init_inst Start the seach from this instance.
     * @param iscorer   the Scoring function.
     * @param max_evals The maximum number of evaluations to perform.
     */
    void operator()(deme_t& deme,
                    const instance& init_inst,
                    const iscorer_base& iscorer,
                    unsigned max_evals,
                    time_t max_time);

    // Like above but assumes that init_inst is null (backward compatibility)
    // In fact, all of the current code uses this entry point, no one
    // bothers to supply an initial instance.
    // XXX This could help PSO if we maintain the best particle.
    void operator()(deme_t& deme,
                    const iscorer_base& iscorer,
                    unsigned max_evals,
                    time_t max_time)
    {
        instance init_inst(deme.fields().packed_width());
        operator()(deme, init_inst, iscorer, max_evals, max_time);
    }

};


} // ~namespace moses
} // ~namespace opencog

#endif
