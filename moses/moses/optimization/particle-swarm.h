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
    ps_parameters(unsigned max_parts = 20,
                // Bit parameters
                double bit_c1 = 0.7,
                double bit_c2 = 1.43,
                // Disc parameters
                double disc_c1 = 2.05,
                double disc_c2 = 2.05,
                // Contin parameters
                double cont_c1 = 0.7,
                double cont_c2 = 1.43
                )
        : max_parts(max_parts),
          bit_c1(bit_c1),
          bit_c2(bit_c2),
          disc_c1(disc_c1),
          disc_c2(disc_c2),
          cont_c1(cont_c1),
          cont_c2(cont_c2),
          inertia_min(0.4),
          inertia_max(0.9),
          bit_min_value(0),
          bit_max_value(1),
          disc_min_value(0),
          disc_max_value(1),
          cont_min_value(std::numeric_limits<contin_t>::min()),
          cont_max_value(std::numeric_limits<contin_t>::max()),
          bit_min_vel(-6), // 1.8%
          bit_max_vel(6), // 98.2%
          disc_min_vel(-disc_max_value/2),
          disc_max_vel(disc_max_value/2),
          cont_min_vel(-std::numeric_limits<contin_t>::max() / 2),
          cont_max_vel(std::numeric_limits<contin_t>::max() / 2),
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

    // For Inertia, i don't know yet what i'll do, but for now i'll maintain
    // one inertia for all types, with values from here:
    // A.P. Engelbrecht, Computational Intelligence: An Introduction, John Willey &
    // Sons Editions, Chichester, 2007.
    // XXX The values probably are those, but the increment factor don't work the
    // same way, it will return some demes before the end of the evaluations.
    double inertia_min, inertia_max; // Min and max inertia weight.

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
    particle_swarm(const optim_parameters& op = optim_parameters())
        : optimizer_base(op), _total_RAM_bytes(getTotalRAM()) {}

protected:
    // Variables:
    const uint64_t _total_RAM_bytes;
    size_t _instance_bytes;
    const ps_parameters ps_params;

    // Functions:
    // log legend for graph stats
    void log_stats_legend();

    void create_random_particle(const field_set& fs,
            instance& new_inst, velocity& vel);

    // Check the limits of something
    void check_bounds(double &val, const double& max, const double& min) {
        if(val < min)
            val = min;
        else if(val > max)
            val = max;
    }

////// Velocity Functions //////
    //// Check bounds functions:
    // There's no real bounds check for bit velocity, it's the probability.
    void check_bit_vel(double &vel) { // Check bounds of bit velocity
        check_bounds(vel, ps_params.bit_min_vel, ps_params.bit_max_vel); }
    void check_disc_vel(double &vel) { // Check bounds of a discrete velocity
        check_bounds(vel, ps_params.disc_min_vel, ps_params.disc_max_vel); }
    void check_cont_vel(double &vel) { // Check bounds of a continuous velocity
        check_bounds(vel, ps_params.cont_min_vel, ps_params.cont_max_vel); }

    //// Generate initial random velocity
    double gen_bit_vel() { // [0,1]
        return (randGen().randdouble() *
                ps_params.range_bit_vel) - ps_params.bit_max_vel; }
    double gen_disc_vel() { // [-0.5, 0.5] when mapped
        return (randGen().randdouble() - ps_params.disc_max_vel); }
    double gen_cont_vel() {
        return (randGen().randdouble() *
                ps_params.range_cont_vel ) - ps_params.cont_max_vel; }

////// Particle values functions //////
    //// Generate initial random instance knob value
    double gen_bit_value() { // [0,1]
        return randGen().randdouble(); }
    double gen_disc_value() { // [0,1]
        return randGen().randdouble(); }
    double gen_cont_value() { // [min, max] of contin_t
        return (((randGen().randbool()) ? 1 : -1) *
                randGen().randdouble() * ps_params.cont_max_value); }

    //// Confinament functions:
    // There isn't confinement for bit values
    // The update rule already kind of do it.
    void confinement_disc(double &value) {
        check_bounds(value, ps_params.disc_min_value, ps_params.disc_max_value); }
    void confinement_cont(double &value) {
        // XXX this will not work, check overflow and underflow
        // before confinement.
        check_bounds(value, ps_params.cont_min_value, ps_params.cont_max_value); }

////// Update specific functions //////

    ////// All
    void update_particle(instance& temp, const instance& local,
            const instance& global, velocity& vels, const field_set& fs);
    //// Bit
    //
    void update_bit_vel(double& vel, int&& temp,
            int&& local, int&& global) { // Bit type is bool
        // Bool convertion to int: false to 0, true to 1.
        // Bit vel hasn't inertia.
        vel += (ps_params.bit_c1 * randGen().randdouble() * (local - temp)) +
            (ps_params.bit_c2 * randGen().randdouble() * (global - temp));
        check_bit_vel(vel);
    }

    // XXX Explanation
    void update_bit_value(bool& value, const double& vel){
        value = (randGen().randdouble() < // Sigmoid
                (1 / (1 + std::exp(-vel)))); // XXX if slow try f(x) = x / (1 + abs(x)) or tanh(x)
    }

    void update_bit_particle(instance& temp, const instance& local, const instance& global,
            const velocity::iterator velocity, const field_set& fs);
    //// Discrete
    // Update discrete velocity
    void update_disc_vel(double& vel, const double& temp,
            const double& local, const double& global) { // Disc type before conversion is double
        // Disc vel hasn't inertia, but has Constriction Factor
        vel += (ps_params.bit_c1 * randGen().randdouble() * (local - temp)) +
            (ps_params.bit_c2 * randGen().randdouble() * (global - temp));
        vel = vel * ps_params.constriction_disc; // Constriction part
        check_disc_vel(vel);
    }

    disc_t cont2disc(const double& cvalue, const unsigned max_dvalue){
        // The original formula is dvalue = std::round(
        // (((max_dvalue - min_dvalue) * (cvalue - min_cvalue))/(max_cvalue - min_cvalue)) + min_dvalue);
        // But [min_cvalue, max_cvalue] == [0,1] and
        // min_dvalue == 0 and max_dvalue == it.multy(), so...
        return (disc_t) std::round(cvalue * max_dvalue); // Return dvalue
    }

    // Discrete values have to be outside the deme.
    struct discrete_particles {
        std::vector<std::vector<double>> temp, best_local;
        unsigned global_index;
        discrete_particles(unsigned part_size, unsigned disc_size) :
            temp(part_size, std::vector<double>(disc_size)),
            best_local(part_size, std::vector<double>(disc_size)),
            global_index(0) {}
    };

    // Update discrete part of the particle
    void update_disc_particle(instance& temp, const instance& local, const instance& global,
            const velocity::iterator velocity, const field_set& fs);

    //// Continuous
    // Update contin part of the particle
    void update_cont_particle(instance& temp, const instance& local, const instance& global,
            const velocity::iterator velocity, const field_set& fs);

    // TODO: Wind dispersion, but test without first
    // Make it later is easy.

////// XXX Remove when new update rules are done.
    struct check_vel_bounds {
        // Create bounds
        // bit: [0,1], disc: [?,?], cont: [-max/2,max/2]
        // XXX hardcoded for now.
    protected:
        void check_bounds(double &val, double max, double min) {
            if(val < min)
                val = min;
            else if(val > max)
                val = max;
        }

    public:
        void bit(double &val) { check_bounds(val, 0, 1);  } // Bit: [0,1]
        void disc(double &val) { check_bounds(val, -4, 4);  } // XXX Disc: [?,?] Probably [-it.multy(), it.multy()]
        void cont(double &val) { check_bounds(val, -100, 100);  } // XXX Cont: [-max/2,max/2]

        double gen_vbit() { return randGen().randdouble(); } // Generate between [0,1]
        double gen_vdisc() { return (randGen().randdouble() * 8) - 4; } // XXX
        double gen_vcont() { return (randGen().randdouble() * 200 ) - 100; } // XXX Cont: [-max/2,max/2]]() {   }; // XXX Cont: [-max/2,max/2]
    } _vbounds;

    // XXX PSO parameters hardcoded just for now.
    //int swarm_size = 20; // Number of particles.
    double cogconst = 0.7, // c1 = Individual learning rate.
    socialconst = 1.43, // c2 = Social parameter.
    inertia_min = 0.4, // wmin = Min of inertia weight.
    inertia_max = 0.9; // wmax = Max of inertia weight.

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
    // XXX In fact, all of the current code uses this entry point, no one
    // bothers to supply an initial instance.
    // This could help PSO if we maintain the best particle.
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
