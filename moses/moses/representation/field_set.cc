/*
 * moses/moses/eda/field_set.cc
 *
 * Copyright (C) 2002-2008 Novamente LLC
 * All Rights Reserved
 *
 * Written by Moshe Looks
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
#include "field_set.h"
#include <sstream>
#include <opencog/util/dorepeat.h>
#include <opencog/util/exceptions.h>
#include <opencog/util/oc_assert.h>
#include <opencog/util/iostreamContainer.h>

namespace opencog {
namespace moses {

const disc_t field_set::term_spec::Stop = 0;

field_set& field_set::operator=(const field_set& rhs)
{
    _fields = rhs._fields;
    _disc = rhs._disc;
    _contin = rhs._contin;
    _term = rhs._term;
    _nbool = rhs._nbool;

    compute_starts();
    return *this;
}
bool field_set::operator==(const field_set& rhs) const
{
    return (_disc == rhs._disc && _contin == rhs._contin &&
            _term == rhs._term && _nbool == rhs._nbool);
}

const term_t& field_set::get_term(const packed_vec& inst, size_t idx) const
{
    size_t raw_idx = term_to_raw_idx(idx);

    // Walk down the tree to get the appropriate term_t.
    const term_spec& o = _term[idx];
    term_tree::iterator it = o.tr->begin();
    for (size_t i = 0; i < o.depth; ++i) {
        disc_t raw_value = get_raw(inst, raw_idx + i);
        if (raw_value == term_spec::Stop) {
            break;
        } else {
            it = o.tr->child(it, term_spec::to_child_idx(raw_value));
        }
    }
    return *it;
}

std::string field_set::to_string(const instance& inst) const
{
    std::stringstream ss;
    ss << "[";
    ostreamContainer(ss, begin_term(inst), end_term(inst), "#", "#", "", false);
    ostreamContainer(ss, begin_contin(inst), end_contin(inst), "|", "|", "", false);
    ostreamContainer(ss, begin_disc(inst), end_disc(inst), " ", " ", "", false);
    ostreamContainer(ss, begin_bit(inst), end_bit(inst), "", "", "", false);
    ss << "]";
    return ss.str();
}

std::string field_set::to_string_raw(const instance& inst) const
{
    std::stringstream ss;
    ostreamContainer(ss, begin_raw(inst), end_raw(inst), "", "[", "]");
    return ss.str();
}

void field_set::build_spec(const spec& s, size_t n)
{
    if (const term_spec* os = boost::get<term_spec>(&s)) {
        build_term_spec(*os, n);
    } else if (const contin_spec* cs = boost::get<contin_spec>(&s)) {
        build_contin_spec(*cs, n);
    } else if (const disc_spec* ds = boost::get<disc_spec>(&s)) {
        build_disc_spec(*ds, n);
    } else {
        OC_ASSERT(false, "This spec is NULL or unknown");
    }
}

void field_set::build_disc_spec(const disc_spec& ds, size_t n)
{
    width_t width = nbits_to_pack(ds.multy);
    size_t base = back_offset();
    for (size_t idx = 0;idx < n;++idx)
        _fields.push_back(field(width,
                                (base + idx*width) / bits_per_packed_t,
                                (base + idx*width) % bits_per_packed_t));
    _disc.insert(_disc.end(), n, ds);
    if (width == 1)
        _nbool += n;
}

void field_set::build_contin_spec(const contin_spec& cs, size_t n)
{
    // All raw contin fields have a multiplicity of 3 (left, right, or
    // stop) and hence are 2 bits wide.  One 'cooked' contin field
    // consists of contin_spec::depth raw fields.
    _contin.insert(_contin.end(), n, cs);
}

void field_set::build_term_spec(const term_spec& os, size_t n)
{
    size_t base = back_offset(), width = nbits_to_pack(os.branching);
    size_t total_width = size_t((width * os.depth - 1) /
                                bits_per_packed_t + 1) * bits_per_packed_t;

    dorepeat (n) {
        dorepeat (os.depth) {
            _fields.push_back(field(width, base / bits_per_packed_t,
                                    base % bits_per_packed_t));
            base += width;
        }
        base += total_width - (os.depth * width); //term vars must pack evenly
    }
    _term.insert(_term.end(), n, os);
}

std::ostream& field_set::ostream_field_set(std::ostream& out) const
{
    using std::endl;

    // Use a pseudo-json style printout.
    out << "field_set = {" << endl;

    out << "n_term_specs= " << term().size()
        << "; n_term_fields= " << n_term_fields()
        << ";\nn_contin_specs= " << contin().size()
        << "; n_contin_fields= " << n_contin_fields()
        << ";\nn_disc_fields= " << n_disc_fields()
        << ";\nn_bit_fields= " << n_bits()
        << ";" << endl;

    unsigned idx = 0;
    out << "fields = {" << endl;

    std::vector<term_spec>::const_iterator tit = term().begin();
    std::vector<term_spec>::const_iterator tend = term().end();
    for (; tit != tend; ++tit, ++idx)
    {
        out << "\t{ idx=" << idx
            << "; type=term"
            << "; depth=" << tit->depth
            << "; branching=" << tit->branching
            << "; }," << endl;
    }

    std::vector<contin_spec>::const_iterator cit = contin().begin();
    std::vector<contin_spec>::const_iterator cend = contin().end();
    for (; cit != cend; ++cit, ++idx)
    {
        out << "\t{ idx=" << idx
            << "; type=contin"
            << "; expansion=" << cit->_exp
            << "; }," << endl;
    }

    std::vector<disc_spec>::const_iterator dit = disc_and_bit().begin();
    std::vector<disc_spec>::const_iterator dend = disc_and_bit().end();
    for (; dit != dend; ++dit, ++idx)
    {
        if (2 < dit->multy)
        {
            out << "\t{ idx=" << idx
                << "; type=disc"
                << "; multiplicity=" << dit->multy
                << "; }," << endl;
        }
        else
        {
            out << "\t{ idx=" << idx
                << "; type=bit; }," << endl;
        }
    }
    out << "}; };" << endl;
    return out;
}


} // ~namespace moses
} // ~namespace opencog
