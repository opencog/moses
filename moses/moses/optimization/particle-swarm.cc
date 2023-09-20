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
#include "../representation/field_set.h"

#include "particle-swarm.h"
#include "hill-climbing.h"

namespace opencog { namespace moses {

////////////////////
// Particle Swarm //
////////////////////

void particle_swarm::operator()(deme_t& best_parts,
                               const instance& init_inst,
                               const iscorer_base& iscorer,
                               unsigned max_evals,
                               time_t max_time)
{
    using namespace boost::placeholders;

    logger().debug("PSO...");

    log_stats_legend();

    // Maintain same name of variables from hill climbing
    // for better understanding.
    // Collect statistics about the run, in struct optim_stats
    nsteps = 0;
    demeID = best_parts.getID();
    over_budget = false;
    struct timeval start;
    gettimeofday(&start, NULL);

    const field_set& fields = best_parts.fields();
    // Track RAM usage. Instances can chew up boat-loads of RAM.
    _instance_bytes = sizeof(instance)
        + sizeof(packed_t) * fields.packed_width();

    unsigned swarm_size = calc_swarm_size(fields);
    // unsigned dim_size = fields.dim_size();

        //if(swarm_size < ps_params.max_parts){
        // If small enough, try all combinations
        // It won't work for if it has contin knobs.
        //    brute_force(best_parts, init_inst, fields, iscorer);
        //    return;
        //}

////// Particle Inicialization //////
    // Reserve uninitialized instances to not have to reallocate.
    // best_parts deme will be the best "personal" (or "local")
    // particles vector to be returned.
    best_parts.reserve(swarm_size);
    // New uninitialized velocity matrix
    // To update the instances
    std::vector<velocity> velocities(swarm_size,
            std::vector<double>(fields.dim_size()));
    // Discrete values of the instance aren't used for update,
    // because of that we need a structure similar to the continuous.
    discrete_particles
        disc_parts(swarm_size, fields.n_disc_fields());
    // Get the 3 uninitialized sets above and initialize then with
    // the instance definition inside fields and ps_params.
    initialize_particles(swarm_size,
            best_parts, velocities, disc_parts, fields);
    // Inicialization of particle best and global best, and their scores.
    auto temp_parts = best_parts;
    unsigned best_global = 0; // Any value
    // Copy the discrete information too
    disc_parts.temp = disc_parts.best_personal;

    // Equal to HC.
    score_t best_score = very_worst_score;
    score_t best_raw_score = very_worst_score;
    size_t current_number_of_evals = 0;

    unsigned iteration = 0;
    unsigned not_improving = 0;
    while(true){
        logger().debug("Iteration: %u", ++iteration);

        // score all instances in the deme
        OMP_ALGO::transform(temp_parts.begin(), temp_parts.end(), temp_parts.begin_scores(),
                            // using bind cref so that score is passed by
                            // ref instead of by copy
                            boost::bind(boost::cref(iscorer), _1));
        current_number_of_evals += swarm_size;

        // XXX What score do i use?
        // I'll use best_score for now.
        bool has_improved = false;
        for (unsigned i = 0; i < swarm_size; ++i) {
            const composite_score &inst_cscore = temp_parts[i].second;
            score_t iscore = inst_cscore.get_penalized_score();
            if(iscore > best_parts[i].second.get_penalized_score()){
                best_parts[i] = temp_parts[i];
                disc_parts.best_personal[i] = disc_parts.temp[i]; //For discrete
                has_improved = true;
                if (iscore > best_score) {
                    best_score = iscore;
                    best_global = i;
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

            // Use raw first than penalized score for find bests
            //const composite_score &inst_cscore = temp_parts[i].second;
            //score_t rscore = inst_cscore.get_score();
            //score_t bp_rscore = best_parts[i].second.get_score();
            //// Comparison
            //if(!(rscore < bp_rscore)){ // Skip if bad
            //    // Get penalizes scores
            //    score_t iscore = inst_cscore.get_penalized_score();
            //    score_t bp_iscore = best_parts[i].second.get_penalized_score();
            //    // Compare
            //    if(rscore != bp_rscore || iscore > bp_iscore){
            //        best_parts[i] = temp_parts[i];
            //        disc_parts.best_personal[i] = disc_parts.temp[i]; //For discrete
            //        has_improved = true;
            //        if (rscore > best_raw_score) {
            //            best_raw_score = rscore;
            //            best_global = i;
            //        }
            //    }
            //}
        }

        // Collect statistics about the run, in struct optim_stats
        struct timeval stop, elapsed;
        gettimeofday(&stop, NULL);
        timersub(&stop, &start, &elapsed);
        start = stop;
        // unsigned usec = 1000000 * elapsed.tv_sec + elapsed.tv_usec;

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
        if (not_improving > 3) {
            logger().debug("Terminate Local Search: Convergence.");
            break;
        }

        // Update particles
        update_particles(temp_parts, best_parts,
                best_global, velocities, disc_parts, fields);
    }

    best_parts.n_best_evals = swarm_size;
    best_parts.n_evals = current_number_of_evals;

} // ~operator

////// The functions below are ordered by utilization order inside the function above.

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

// TODO: Explanation
// There's no explanation for this, it's just a temporary solution.
// Maybe use adaptative pso, something like LPSO (Lander).
unsigned particle_swarm::calc_swarm_size(const field_set& fs) {
    // For disc i'll use the same proportion of bit.
    const double byte_relation = 3.0 / 4.0, // For each 4 bit i'll let it similar to a cont.
                cont_relation = 3; // Normally 3x or 4x of the dimension.

    unsigned disc_bit_size = fs.n_bits();

    double total = disc_bit_size * byte_relation +
                    fs.contin().size() * cont_relation;
    total = clamp(total, 4.0, (double) ps_params.max_parts); // 4 For min, less than this is almost useless.
    return std::round(total); // Round it.
}

//void particle_swarm::brute_force(deme_t deme, const instance& init_inst,
//        const field_set& fs, const iscorer_base& iscorer){
//    // Add init_inst
//    deme.push_back(init_inst);
//    // Generate all possibilities
//    generate_all_in_neighborhood(fs, fs.dim_size(),
//            deme.begin() + 1, //Plus the init_inst
//            deme.end(), init_inst);
//
//    // score all instances in the deme
//    OMP_ALGO::transform(deme.begin(), deme.end(),
//            deme.begin_scores(), boost::bind(boost::cref(iscorer), _1));
//
//    deme.n_best_evals = deme.size();
//    deme.n_evals = deme.size();
//}

void particle_swarm::initialize_particles (const unsigned& swarm_size,
        deme_t& best_parts, std::vector<velocity>& velocities,
        discrete_particles& disc_parts, const field_set& fields) {
    auto velit = velocities.begin();
    auto dvit = disc_parts.best_personal.begin();
    dorepeat(swarm_size){
        instance new_inst(fields.packed_width());
        initialize_random_particle(
            new_inst, *velit, *dvit, fields);
        velit++; dvit++;
        best_parts.push_back(new_inst);
    }
}

void particle_swarm::initialize_random_particle (instance& new_inst,
        velocity& vel, std::vector<double>& dist_values, const field_set& fs){
    auto vit = vel.begin();
    // For each bit
    for(auto it = fs.begin_bit(new_inst);
            it != fs.end_bit(new_inst); ++it, ++vit) {
        *it = gen_bit_value(); // New bit value in instance
        *vit = gen_bit_vel(); // New bit velocity
    }
    // For each disc
    auto cit = dist_values.begin();
    for(auto dit = fs.begin_disc(new_inst);
            dit != fs.end_disc(new_inst); ++dit, ++vit, ++cit) {
        *cit = gen_disc_value(); // New cont value for disc
        *dit = cont2disc(*cit, dit.multy()); // New disc value in instance
        *vit = gen_disc_vel(); // New disc velocity
    }
    // For each contin
    for(auto it = fs.begin_contin(new_inst);
            it != fs.end_contin(new_inst); ++it, ++vit) {
        *it = gen_cont_value(); // New cont value in instance
        *vit = gen_cont_vel(); // New cont velocity
    }
}

void particle_swarm::update_particles(deme_t& temp_parts, const deme_t& best_parts, const int& best_index,
        std::vector<velocity>& velocities, discrete_particles& disc_parts, const field_set& fields) {

    // Iteration initialization
    auto temp_it = temp_parts.begin();
    auto end_temp = temp_parts.end(); // End
    auto bestp_it = best_parts.begin();
    auto vels_it = velocities.begin();
    // Discrete part
    auto disc_temp_it = disc_parts.temp.begin();
    auto disc_best_it = disc_parts.best_personal.begin();

    // Get best instance from best particles
    const instance& best_instance = best_parts[best_index].first;
    const std::vector<double>& disc_best_instance =
                            disc_parts.best_personal[best_index];

    // Over each instance
    for(;temp_it != end_temp; temp_it++, bestp_it++,
            vels_it++, disc_temp_it++){
        instance& temp_inst = (*temp_it).first;
        const instance& bestp_inst = (*bestp_it).first;
        velocity::iterator vel = (*vels_it).begin();

        update_bit_particle(temp_inst, bestp_inst, best_instance, vel, fields);
        update_disc_particle(temp_inst, *disc_temp_it,
                *disc_best_it, disc_best_instance, vel, fields);
        update_cont_particle(temp_inst, bestp_inst, best_instance, vel, fields);
    }
}

void particle_swarm::update_bit_particle(instance& temp, const instance& personal,
        const instance& global, velocity::iterator vel, const field_set& fs) {
    auto tit = fs.begin_bit(temp); // Iterator for temporary particle
    auto pit = fs.begin_bit(personal), // Iterator for best personal particle
         git = fs.begin_bit(global); // Iterator for best instance
    for(;tit != fs.end_bit(temp); // Comparison
            ++tit, ++pit, ++git, ++vel){ //Next
        update_bit_vel(*vel, *tit, *pit, *git);
        *tit = new_bit_value(*vel);
    }
}

void particle_swarm::update_disc_particle(instance& dtemp, std::vector<double>& temp,
        const std::vector<double>& personal, const std::vector<double>& global,
        velocity::iterator vel, const field_set& fs) {
    auto disc_it = fs.begin_disc(dtemp); // Iterator for temporary discrete particle
    auto tit = temp.begin(); // Iterator for temporary particle
    auto pit = personal.begin(), // Iterator for best personal particle
         git = global.begin(); // Iterator for best instance
    for(;disc_it != fs.end_disc(dtemp); // Comparison
            disc_it++, ++tit, ++pit, ++git, ++vel){ //Next
        update_disc_vel(*vel, *tit, *pit, *git);
        *disc_it = new_disc_value(*tit, *vel, disc_it.multy());
    }
}

void particle_swarm::update_cont_particle(instance& temp, const instance& personal,
        const instance& global, velocity::iterator vel, const field_set& fs) {
    auto tit = fs.begin_contin(temp); // Iterator for temporary particle
    auto pit = fs.begin_contin(personal), // Iterator for best personal particle
         git = fs.begin_contin(global); // Iterator for best instance
    for(;tit != fs.end_contin(temp); // Comparison
            ++tit, ++pit, ++git, ++vel){ //Next
        update_cont_vel(*vel, *tit, *pit, *git);
        *tit = new_cont_value(*tit, *vel);
    }
}


} // ~namespace moses
} // ~namespace opencog

