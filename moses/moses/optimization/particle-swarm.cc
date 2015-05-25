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
    logger().debug("Enter PSO...");
// TODO: basically everything.
    optimizer_base* _opthc =
        new hill_climbing(opt_params, ps_params);
    optimizer_base& _hc = *_opthc;
    _hc(deme, iscorer, max_evals, max_time);
#define ACCEPTABLE_SIZE 5000
#define ACCEPTABLE_RAM_FRACTION 0.5

} // ~operator

// TODO: verify comments and code
/// Shrink the deme, by removing all instances with score less than
/// 'cutoff'.  This is implemented with in-place deletion of elements
/// from a vector, with at least a token attempt to delete contigous
/// regions of low scores, in one go.  It is possible that a faster
/// algorithm would be to sort first, and then delete the tail-end of
/// the vector.  But this ixsn't know ... XXX experiment with this!?
/// ... err, but right now, trimming takes a small fraction of a second,
/// so there is no rush to fis this.
size_t particle_swarm::resize_by_score(deme_t& deme, score_t cutoff)
{
    size_t ndeleted = 0;
    while (true) {
         auto first = deme.end();
         auto last = deme.end();
         size_t contig = 0;
         for (auto it = deme.begin(); it != deme.end(); it++) {
             score_t iscore = it->second.get_penalized_score();
             if (iscore <= cutoff) {
                 if (0 == contig) first = it;
                 last = it;
                 contig ++;
             } else {
                 if (0 < contig)
                     break;
             }
         }

         if (0 == contig)
             break;

         if (last != deme.end()) last++;
         deme.erase(first, last);
         ndeleted += contig;

         // Keep around at least 20 instances; useful for cross-over.
#define MIN_TO_KEEP 20
         if (deme.size() < MIN_TO_KEEP)
             break;
    }
    logger().debug() << "Trimmed "
            << ndeleted << " low scoring instances.";
    return deme.size();
}

/// Keep the size of the deme at a managable level.
/// Large populations can easily blow out the RAM on a machine,
/// so we want to keep it at some reasonably trim level.
//
bool particle_swarm::resize_deme(deme_t& deme, score_t best_score)
{
    bool did_resize =  false;
    // Lets see how many we might be able to trounce.
    score_t cutoff = best_score - ps_params.score_range;
    size_t bad_score_cnt = 0;

    // To find the number of bad scores, we have to look
    // at the *whole* deme.
    for (const deme_inst_t& si : deme) {
        score_t iscore = si.second.get_penalized_score();
        if (iscore <=  cutoff)
            bad_score_cnt++;
    }


    // To avoid wasting cpu time pointlessly, don't bother with
    // population size management if we don't get any bang out
    // of it.
#define DONT_BOTHER_SIZE 500
    uint64_t usage = _instance_bytes * deme.size();

    if ((DONT_BOTHER_SIZE < bad_score_cnt) or
        (ACCEPTABLE_RAM_FRACTION * _total_RAM_bytes < usage))
    {

        logger().debug() << "Will trim " << bad_score_cnt
            << " low scoring instances out of " << deme.size();
        resize_by_score(deme, cutoff);
        did_resize = true;
    }

    // Are we still too large? Whack more, if needed.
    // We want to whack only the worst scorerers, and thus
    // a partial sort up front.
    if ((ps_params.max_allowed_instances < deme.size()) or
        (ACCEPTABLE_RAM_FRACTION * _total_RAM_bytes < usage))
    {
        std::partial_sort(deme.begin(),
                          next(deme.begin(), ps_params.max_allowed_instances),
                          deme.end(),
                          std::greater<deme_inst_t>());

        deme.erase(next(deme.begin(), ps_params.max_allowed_instances),
                   deme.end());
        did_resize = true;
    }
    return did_resize;
}

void particle_swarm::log_stats_legend()
{
    logger().info() << ps_params.prefix_stat_deme << "Hill: # "   /* Legend for graph stats */
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

