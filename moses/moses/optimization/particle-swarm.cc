/*
 * moses/moses/optimization/particle-swarm.cc
 *
 * Copyright (C) 2002-2008 Novamente LLC
 * Copyright (C) 2012 Poulin Holdings LLC
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

#include <math.h>   // for sqrtf, cbrtf

#include <boost/algorithm/minmax_element.hpp>

#include <opencog/util/oc_omp.h>

#include "../moses/neighborhood_sampling.h"

#include "particle-swarm.h"

namespace opencog { namespace moses {

////////////////////
// Particle Swarm //
////////////////////

void particle_swarm::operator()(deme_t& deme,
                               const instance& init_inst,
                               const iscorer_base& iscorer,
                               unsigned max_evals,
                               time_t max_time)
{
    logger().debug("PSO...");

    log_stats_legend();

    // Maintain same name of variables from hill climbing
    // for better understanding.
    // Collect statistics about the run, in struct optim_stats
    nsteps = 0;
    demeID = deme.getID();
    over_budget = false;
    struct timeval start;
    gettimeofday(&start, NULL);

    const field_set& fields = deme.fields();
    // Track RAM usage. Instances can chew up boat-loads of RAM.
    _instance_bytes = sizeof(instance)
        + sizeof(packed_t) * fields.packed_width();

    // Swarm size has to be variable, it would be a shame to use a lot of
    // particles when you don't need (dim == 1);
    int swarm_size = 4;
    int dim_size = fields.dim_size();
    // Inertia calculation
    double inertia_factor = inertia_max,
           decinertia_factor = (inertia_max - inertia_min) / (max_evals / swarm_size);

    // TODO: Update velocity if deme.size() < swarm_size
    // I don't know yet how can i know what instances were throw away.
    // inheritance didn't works because it isn't a vector of points.
    //for(auto inst : deme)
    //    update_vel(inst._vel);
    std::vector<velocity> _velocities;
    dorepeat(deme.size()){
        velocity vel(dim_size);
        // Reinicialize new velocities.
        dorepeat(fields.n_bits()) { vel.push_back(_vbounds.gen_vbit()); }
        dorepeat(fields.n_disc_fields()) { vel.push_back(_vbounds.gen_vdisc()); }
        dorepeat(fields.n_contin_fields()) { vel.push_back(_vbounds.gen_vcont()); }
        _velocities.push_back(vel);
    }
    // Deme size == particle size.
    dorepeat(swarm_size - deme.size()){
        instance new_inst(fields.packed_width());
        velocity vel(dim_size);
        create_random_particle(fields, new_inst, vel);
        deme.push_back(new_inst);
        _velocities.push_back(vel);
    }

    // Inicialization of particle best and global best, and their scores.
    auto bests_part = deme;
    auto best_global = deme[0];
    // Equal to HC.
    composite_score best_cscore = worst_composite_score;
    score_t best_score = very_worst_score;
    score_t best_raw_score = very_worst_score;
    size_t current_number_of_evals = 0;

    unsigned iteration = 0;
    unsigned not_improving = 0;
    while(true){
        logger().debug("Iteration: %u", ++iteration);

        // score all instances in the deme
        OMP_ALGO::transform(deme.begin(), deme.end(), deme.begin_scores(),
                            // using bind cref so that score is passed by
                            // ref instead of by copy
                            boost::bind(boost::cref(iscorer), _1));

        // XXX What score do i use?
        // I'll use best_score for now.
        bool has_improved = false;
        for (unsigned i = 0; i < deme.size(); ++i) {
            const composite_score &inst_cscore = deme[i].second;
            score_t iscore = inst_cscore.get_penalized_score();
            if(iscore > bests_part[i].second.get_penalized_score()){
                bests_part[i] = deme[i];
                if (iscore >  best_global.second.get_penalized_score()) {
                    has_improved = true;
                    best_cscore = inst_cscore;
                    best_score = iscore;
                    best_global = deme[i];
                }
            }

            // The instance with the best raw score will typically
            // *not* be the same as the the one with the best
            // weighted score.  We need the raw score for the
            // termination condition, as, in the final answer, we
            // want the best raw score, not the best weighted score.
            score_t rscore = inst_cscore.get_score();
            if (rscore >  best_raw_score) {
                best_raw_score = rscore;
            }
        }

        current_number_of_evals += swarm_size;

        // Collect statistics about the run, in struct optim_stats
        struct timeval stop, elapsed;
        gettimeofday(&stop, NULL);
        timersub(&stop, &start, &elapsed);
        start = stop;
        unsigned usec = 1000000 * elapsed.tv_sec + elapsed.tv_usec;

        /* If we've blown our budget for evaluating the scorer,
         * then we are done. */
        if (max_evals <= current_number_of_evals) {
            over_budget = true;
            logger().debug("Terminate Local Search: Over budget");
            break;
        }

        if (max_time <= elapsed.tv_sec) {
            over_budget = true;
            logger().debug("Terminate Local Search: Out of time");
            break;
        }
        max_time -= elapsed.tv_sec; // count-down to zero.

        /* If we've aleady gotten the best possible score, we are done. */
        if (opt_params.terminate_if_gte <= best_raw_score) {
            logger().debug("Terminate Local Search: Found best score");
            break;
        }

        // TODO: work in a better way to identify convergence.
        not_improving = (has_improved) ? 0 : not_improving + 1;
        if (not_improving > 1) {
            logger().debug("Terminate Local Search: Convergence.");
            break;
        }

        // TODO: update velocities.
        for(unsigned part = 0; part < swarm_size; ++part){ // Part == particle index
            unsigned dim = 0; // Dim == dimension index

            // Bit velocity update
            // XXX IT ISN'T THE ORIGINAL BPSO it's just a test for now
            for(auto dit = fields.begin_bit(deme[part].first),
                    lit = fields.begin_bit(bests_part[part].first),
                    git = fields.begin_bit(best_global);
                    dit != fields.end_bit(deme[part].first);
                    ++dit, ++lit, ++git, ++dim){ //Next

                double& vel = _velocities[part][dim];
                vel = (inertia_factor * vel) + // Maintain velocity
                    (cogconst * randGen().randdouble() * (*lit - *dit)) + // Go to my best
                    (socialconst * randGen().randdouble() * (*git - *dit)); // Go to global best
                _vbounds.bit(vel); // check bounds for bit velocity
                //*dit = ((*dit + vel) > 0.5) ? true : false;
                //logger().debug("test: %d", *dit);
            }
            // Discrete velocity update
            // XXX IT ISN'T THE ORIGINAL DPSO it's just a test for now
            for(auto dit = fields.begin_disc(deme[part].first),
                    lit = fields.begin_disc(bests_part[part].first),
                    git = fields.begin_disc(best_global);
                    dit != fields.end_disc(deme[part].first);
                    ++dit, ++lit, ++git, ++dim){ //Next

                double& vel = _velocities[part][dim];
                vel = (inertia_factor * vel) + // Maintain velocity
                    (cogconst * randGen().randdouble()
                     * (double)(*lit - *dit)) + // Go to my best
                    (socialconst * randGen().randdouble()
                     * (double)(*git - *dit)); // Go to global best
                _vbounds.disc(vel); // check bounds for disc velocity
                *dit = std::fmin(std::fmax((((double) *dit) + vel), 7), 0);// CONFINAMENT WITHOUT WIND DISPERSION
            }
            // Continuos velocity update
            // This is the original one
            for(auto dit = fields.begin_contin(deme[part].first),
                    lit = fields.begin_contin(bests_part[part].first),
                    git = fields.begin_contin(best_global);
                    dit != fields.end_contin(deme[part].first);
                    ++dit, ++lit, ++git, ++dim){ //Next

                double& vel = _velocities[part][dim];
                vel = (inertia_factor * vel) + // Maintain velocity
                    (cogconst * randGen().randdouble() * (*lit - *dit)) + // Go to my best
                    (socialconst * randGen().randdouble() * (*git - *dit)); // Go to global best
                _vbounds.cont(vel); // check bounds for contin velocity
                //*dit = *dit + vel;
            }
        }

        // TODO: update particles
        inertia_factor -= decinertia_factor;
    }

    deme.n_best_evals = swarm_size;
    deme.n_evals = current_number_of_evals;

} // ~operator

void particle_swarm::create_random_particle (
        const field_set& fs, instance& new_inst, velocity& vel){
    // For each bit
    if(fs.n_bits() > 0) {
        for(auto it = fs.begin_bit(new_inst);
                it != fs.end_bit(new_inst); ++it) {
            *it = randGen().randbool();
            vel.push_back(_vbounds.gen_vbit());
        }
    }
    // For each disc
    if(fs.n_disc_fields() > 0) {
        for(auto it = fs.begin_disc(new_inst);
                it != fs.end_disc(new_inst); ++it) {
            *it = randGen().randint(it.multy());
            vel.push_back(_vbounds.gen_vdisc());
        }
    }
    // For each bit
    if(fs.n_contin_fields() > 0) {
        for(auto it = fs.begin_contin(new_inst);
                it != fs.end_contin(new_inst); ++it) {
            *it = randGen().randint();
            vel.push_back(_vbounds.gen_vcont());
        }
    }
}

void particle_swarm::log_stats_legend()
{
    logger().info() << "PSO: # "   /* Legend for graph stats */
        "demeID\t"
        "iteration\t"
        "total_steps\t"
        "total_evals\t"
        "microseconds\t"
        "new_instances\t"
        "num_instances\t"
        "inst_RAM\t"
        "num_evals\t"
        "has_improved\t"
        "best_weighted_score\t"
        "delta_weighted\t"
        "best_raw\t"
        "delta_raw\t"
        "complexity";
}

} // ~namespace moses
} // ~namespace opencog

