/** interpreter.h --- 
 *
 * Copyright (C) 2012 OpenCog Foundation
 *
 * Author: Nil Geisweiller
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


#ifndef _OPENCOG_INTERPRETER_H
#define _OPENCOG_INTERPRETER_H

#include "../combo/vertex.h"

/**
 * Attempt to rewrite the combo interpreter into something modular and efficient.
 *
 * There are several interpreter each dealing with its data type. One
 * can use a sub-interpreter for a specific data type at a lower cost
 * of storing the variable content and avoiding the overhead of
 * building a vertex or a combo_tree.
 */

namespace opencog { namespace combo {

/**
 * Interpreter for boolean expressions.
 *
 * We use builtin as inputs and output because boolean contants are
 * encoded as builtins id::logical_true id::logical_false.
 */
struct boolean_interpreter
{
    // ctor
    boolean_interpreter(const std::vector<builtin>& inputs = std::vector<builtin>());

    // interpreters
    builtin operator()(const combo_tree& tr) const;
    builtin operator()(const combo_tree::iterator it) const;
    virtual builtin boolean_eval(combo_tree::iterator it) const;

protected:
    const std::vector<builtin>& boolean_inputs;
};

/**
 * Interpreter for contin expressions.
 *
 * We use contin as inputs and output
 */
struct contin_interpreter
{
    // ctor
    contin_interpreter(const std::vector<contin_t>& inputs = std::vector<contin_t>());

    // interpreters
    contin_t operator()(const combo_tree& tr) const;
    contin_t operator()(const combo_tree::iterator it) const;
    virtual contin_t contin_eval(combo_tree::iterator it) const;

protected:
    const std::vector<contin_t>& contin_inputs;
};

struct mixed_interpreter : public boolean_interpreter, public contin_interpreter
{
    // ctor
    mixed_interpreter(const std::vector<vertex>& inputs = std::vector<vertex>());
    mixed_interpreter(const std::vector<contin_t>& inputs);
    mixed_interpreter(const std::vector<builtin>& inputs);

    // interpreters
    vertex operator()(const combo_tree& tr) const;
    vertex operator()(const combo_tree::iterator it) const;
    virtual builtin boolean_eval(combo_tree::iterator it) const;
    virtual contin_t contin_eval(combo_tree::iterator it) const;
    virtual vertex mixed_eval(combo_tree::iterator it) const;

protected:
    bool _use_boolean_inputs;
    bool _use_contin_inputs;
    const std::vector<vertex>& _mixed_inputs;
};            

// TODO: port the functional interpreter to support list and lambda
#if 0
    if (const builtin* b = boost::get<builtin>(&v)) {
        switch (*b) {

        // list constructor
        case id::list : {

            combo_tree tr(id::list);
            pre_it loc = tr.begin();

            for (sib_it sib = it.begin(); sib != it.end(); ++sib) {
                // tr.append_child(loc, eval_throws_tree(bmap, sib).begin());
                combo_tree rr = eval_throws_tree(bmap, sib);
                tr.append_child(loc, rr.begin());
             }

            return tr;
        }

        // car takes a list and returns head of the list
        case id::car : {
            sib_it lp = it.begin();

            combo_tree evo;
            if (*lp != id::list) {
                evo = eval_throws_tree(bmap, lp);
                lp = evo.begin();
            }
            if (*lp != id::list)
                throw ComboException(TRACE_INFO, "not a list!");

            // If the list is empty; then return empty list!
            // That is, use an empty list to represent nil.
            if (lp.begin() == lp.end())
                return combo_tree(id::list);
            return eval_throws_tree(bmap, lp.begin());
        }

        // cdr takes a list and returns the tail of the list
        case id::cdr : {
            sib_it top = it.begin();

            combo_tree evo;
            if (*top != id::list) {
                evo = eval_throws_tree(bmap, top);
                top = evo.begin();
            }
            if (*top != id::list)
                throw ComboException(TRACE_INFO, "not a list!");

            sib_it sib = top.begin();

            // Skip over the first elt
            sib ++;

            combo_tree tr(id::list);
            pre_it loc = tr.begin();
            for (; sib != top.end(); ++sib) {
                combo_tree rest = eval_throws_tree(bmap, sib);
                tr.append_child(loc, rest.begin());
            }

            return tr;
        }

        // cons takes an element and a list,
        // and adds the element to the head of the list
        case id::cons : {
            combo_tree tr(id::list);
            pre_it loc = tr.begin();
            if(it.begin()==it.end()){
              vertex vx(id::cons);
              return combo_tree(vx);
            }
            sib_it head = it.begin();
            combo_tree ht = eval_throws_tree(bmap, head);
            tr.append_child(loc, ht.begin());

            ++head;
            combo_tree rest = eval_throws_tree(bmap, head);

            sib_it lst = rest.begin();
            for (sib_it sib = lst.begin(); sib != lst.end(); ++sib)
                // tr.append_child(loc, eval_throws_tree(bmap, sib).begin());
                tr.append_child(loc, (pre_it) sib);

            return tr;
        }

        case id::foldr : {
            combo_tree tr(it);
 
            // base case: foldr(f v list) = v
            // i.e. list is empty list.
            sib_it itend = tr.begin().end();
            --itend;
            if (itend.begin() == itend.end()) {
                --itend;
                return eval_throws_tree(bmap, itend);
            }

            // main case: foldr(f v l) = f(car foldr(f v cdr))
            // new tree: f(car foldr(f v cdr))

            // cb_tr will be the call to f.
            sib_it f = it.begin();
            combo_tree cb_tr(f);
            sib_it loc = cb_tr.begin();

            // eval: car(<the list>)
            combo_tree car_lst(id::car); 
            sib_it car_lst_it = car_lst.begin();
            car_lst.append_child(car_lst_it, itend);
            car_lst = eval_throws_tree(bmap,car_lst_it);
            car_lst_it = car_lst.begin();

            // cb_tr == f(<x>)
            cb_tr.append_child(loc, car_lst_it);

            // copy {the combo subtree for this use of foldr} and change it to include the cdr of L instead of L
            // let cdr_call = a new tree containing a call to cdr
            combo_tree cdr_call(id::cdr);
            sib_it cdr_call_it = cdr_call.begin();
            cdr_call.append_child(cdr_call_it, itend);
            tr.erase(itend);
            sib_it tr_it = tr.begin();

            // let cdr_result = the result of the call to cdr
            combo_tree cdr_result(eval_throws_tree(bmap, cdr_call_it));
            sib_it cdr_result_it = cdr_result.begin();
            tr.append_child(tr_it, cdr_result_it);

            tr = eval_throws_tree(bmap, tr_it);
            tr_it = tr.begin();

            cb_tr.append_child(loc, tr_it);

            return eval_throws_tree(bmap, loc);
        }

        case id::foldl : {
            combo_tree tr(it);

            // base case: foldl(f v list) = v
            // i.e. list is empty list.
            sib_it itend = tr.begin().end();
            --itend;
            if (itend.begin() == itend.end()) {
                --itend;
                return eval_throws_tree(bmap, itend);
            }

            // new tree: foldl(f f(v car) cdr)
            sib_it f = it.begin();
            combo_tree lst = eval_throws_tree(bmap, itend);
            sib_it lst_it = lst.begin();
            sib_it tr_it = tr.begin();
            tr.erase(itend);
            --itend;
            tr.erase(itend);
	    
            combo_tree rec(f);
            sib_it rec_it = rec.begin();
            sib_it v = ++f;
            rec.append_child(rec_it, v);
            combo_tree car_lst(id::car);
            sib_it car_lst_it = car_lst.begin();
            car_lst.append_child(car_lst_it, lst_it);
            car_lst = eval_throws_tree(bmap, car_lst_it);
            car_lst_it = car_lst.begin();
            rec.append_child(rec_it, car_lst_it);  //f(v car)
            tr.append_child(tr_it, rec_it);

            combo_tree cdr_lst(id::cdr);
            sib_it cdr_lst_it = cdr_lst.begin();
            cdr_lst.append_child(cdr_lst_it, lst_it);
            cdr_lst = eval_throws_tree(bmap, cdr_lst_it);
            cdr_lst_it = cdr_lst.begin();
            tr.append_child(tr_it, cdr_lst_it);

            return eval_throws_tree(bmap, tr_it ); 
        }

        // lambda constructor
        case id::lambda : {
            combo_tree tr(it);
            return tr;
        }

        // The apply() operator is a sensible thing to have, but this code is bad so I disabled it.
        // It shouldn't use set_bindings; if we want lambda functions then we should use scopes properly -- Jade
        case id::apply : {
            OC_ASSERT(false, "apply() is not implemented");
            combo_tree tr(it);
            sib_it tr_it = tr.begin().begin();
            sib_it lambda_it = tr_it.end();
            lambda_it --;
            combo_tree exp_tr(lambda_it);

            vector<vertex> al; // list of arguments
            ++tr_it;
            for(; tr_it!=tr.begin().end(); ++tr_it){
                al.push_back(*tr_it);
            }
            //set_bindings(exp_tr, exp_tr.begin(), al, explicit_arity(exp_tr));
            return eval_throws_tree(bmap, exp_tr);
        }
#endif		
		
}} // ~namespaces combo opencog
        
#endif // _OPENCOG_INTERPRETER_H
