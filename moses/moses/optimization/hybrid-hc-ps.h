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
#ifndef _MOSES_HYBRID_HC_PS_H
#define _MOSES_HYBRID_HC_PS_H

#include <opencog/util/oc_assert.h>

#include "../representation/instance_set.h"
#include "optimization.h"
#include "hill-climbing.h"
#include "particle-swarm.h"

namespace opencog { namespace moses {

typedef std::vector<double> velocity;

struct hybrid_hc_ps : optimizer_base
{
    hybrid_hc_ps(const optim_parameters& op = optim_parameters(),
                    const hc_parameters& hc = hc_parameters(),
                    const ps_parameters& ps = ps_parameters())
        : optimizer_base(op), _total_RAM_bytes(getTotalRAM()), hc_params(hc)
    {
        get_params(ps);
    }

protected:
    struct hybrid_parameters {
        // From Particle Swarm
        double _min_vel, _max_vel, // Max and min velocities (normally min/2,max/2 values)
               _min_value, _max_value, // Max and min for contin values (for now dependent of contin depth)
               _range_vel, _range_value, // Range Values (max - min)
               _inertia, _c1, _c2; // Update parameters
    };

    // Get the necessary parameters from hill climbing and from particle swarm.
    void get_params(const ps_parameters& ps){
        // Particle Swarm
        // I will use 2^depth as max value to match how neighborhood is using it now.
        ps_params._max_value = pow2(get_depth());
        ps_params._min_value = -ps_params._max_value;
        ps_params._max_vel = ps_params._max_vel / 2;
        ps_params._min_vel = ps_params._min_vel / 2;
        ps_params._range_value = ps_params._max_value - ps_params._min_value;
        ps_params._range_vel = ps_params._max_vel - ps_params._min_vel;

        ps_params._inertia = ps.inertia;
        ps_params._c1 = ps.cont_c1;
        ps_params._c2 = ps.cont_c2;
    }

    // Variables:
    const uint64_t _total_RAM_bytes;
    size_t _instance_bytes;
    const hc_parameters hc_params;
    hybrid_parameters ps_params;

//////////////////////////////
// Hill Climbing Functions //
//////////////////////////////

    // log legend for graph stats
    void log_stats_legend();

    // Return an estimate of the size of the neighborhood at distance
    // 'distance'
    size_t estimate_neighborhood(size_t distance, const field_set& fields);

    // Return an estimate of the number of new instances to search
    size_t n_new_instances(size_t distance, unsigned max_evals,
                           size_t current_number_of_evals,
                           size_t total_number_of_neighbors);

    /**
     * Cross the single top-scoring instance against the next-highest scorers.
     *
     * As arguments, accepts a range of scored instances ("the sample"),
     * and a single instance from which these were all derived ("the base"
     * or center instance).  This will create a number of new instances,
     * which will be a cross of the highest-scoring instance with the
     * next-highest scoring instances.
     *
     * @deme:         the deme holding current instances, and where
     *                new instances will be placed.
     * @deme_size:    the current size of the deme. New instances
     *                will be appended at the end.
     * @base:         the base instance from which the sample was
     *                was derived.
     * @sample_start: the count, within the deme, at which the
     *                scored instances start. These are assumed to
     *                have been derived from the base instance.
     * @sample_size:  the number of instances in the sample. These
     *                are assumed to be in sequential order, starting
     *                at sample_start.
     * @num_to_make:  Number of new instances to create.  The actual
     *                number created will be the lesser of this and
     *                sample_size-1.
     */
    size_t cross_top_one(deme_t& deme,
                         size_t deme_size,
                         size_t num_to_make,
                         size_t sample_start,
                         size_t sample_size,
                         const instance& base);

    /** two-dimensional simplex version of above. */
    size_t cross_top_two(deme_t& deme,
                         size_t deme_size,
                         size_t num_to_make,
                         size_t sample_start,
                         size_t sample_size,
                         const instance& base);

    /** three-dimensional simplex version of above. */
    size_t cross_top_three(deme_t& deme,
                           size_t deme_size,
                           size_t num_to_make,
                           size_t sample_start,
                           size_t sample_size,
                           const instance& base);

    // chain the 3 crossovers methods above and return the number of new instances
    size_t crossover(deme_t& deme, size_t deme_size,
                     size_t sample_start, size_t sample_size,
                     const instance& base);

    bool resize_deme(deme_t& deme, score_t score_cutoff);
    size_t resize_by_score(deme_t& deme, score_t score_cutoff);

//////////////////////////////
// Particle Swarm Functions //
//////////////////////////////

////// Velocity Functions //////
    //// Check bounds functions:
    void check_cont_vel(double &vel) { // Check bounds of a continuous velocity
        vel = bound(vel, ps_params._min_vel, ps_params._max_vel); }

    //// Generate initial random velocity
    double gen_cont_vel() {
        return (randGen().randdouble() *
                ps_params._range_vel) - ps_params._max_vel;
    }

////// Particle values functions //////
    //// Generate initial random instance knob value
    double gen_cont_value() { //
        return (randGen().randdouble() *
                ps_params._range_value) + ps_params._min_value;
    }

    //// Confinament functions:
    void confinement_cont(double& value) {
        value = bound(value, ps_params._min_value,
                ps_params._max_value);
    }

////// Update specific functions //////
    //// All
    //void update_particles(deme_t& temp_parts, const deme_t& best_parts, const int& best_index,
    //    std::vector<velocity>& velocities, discrete_particles& disc_parts, const field_set& fields);

    //// Continuous
    // Update continuous velocity
    void update_cont_vel(double& vel, const double& temp,
            const double& personal, const double& global) { // Cont type is double
        vel += ps_params._inertia * ((ps_params._c1 * randGen().randdouble() * (personal - temp)) +
                (ps_params._c2 * randGen().randdouble() * (global - temp)));
        check_cont_vel(vel);
    }

    // XXX Explanation
    contin_t new_cont_value(contin_t& value, const double& vel){
        // Wind dispersion enter here.
        // XXX Check overflow
        value += vel;
        confinement_cont(value);
        return value;
    }

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
        instance init_inst(deme.fields().packed_width(), deme.fields().n_contin_fields());
        //XXX randomize
        operator()(deme, init_inst, iscorer, max_evals, max_time);
    }

};


} // ~namespace moses
} // ~namespace opencog

#endif
