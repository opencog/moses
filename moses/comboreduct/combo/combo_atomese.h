//
// Created by kasim on 5/19/18.
//

#ifndef MOSES_COMBO_ATOMESE_H
#define MOSES_COMBO_ATOMESE_H

#include <opencog/atoms/base/Handle.h>
#include "vertex.h"

namespace opencog{namespace combo{
        using namespace opencog;
        /**
         * Convert a combo_tree to atomese program.
         *
         * @param combo_tree   the combo_tree containing the combo program
         * @return                 the Handle containing the atomese program
         */
        Handle atomese_combo(const combo_tree &tr);

        /**
         * Convert a combo_tree to atomese from a head of a combo_tree program.
         *
         * @param combo_tree::iterator   the iterater to the head of a combo_tree
         * @return                 the Handle containing the atomese program
         */
        template<typename Iter>
        Handle atomese_combo_it(Iter, id::procidure_types&);

        /**
         * Convert a combo_tree::vertex to atomese program.
         *
         * @param Handle&       handle ref, if the vertex is converted to node it will be stored here
         * @param const vertex&     a vertex ref containing the combo_tree vertex to be converted to atom
         * @param procidure_type&     containing the type  of the atom i:e predicate, shema
         * @return                 if the vertex is to be converted to an atomese link return Link type
         *                          otherwise return -1
         */
        Type atomese_vertex(Handle&, const vertex&, id::procidure_types&);

        /**
         * Convert a combo argument to atomese.
         *
         * @param Handle&       handle ref, stores the atomese converted from combo argument
         * @param const argument&     a argumet ref containing the combo argument to be converted to atom
         * @param procidure_type&     ref to parent procidere type containing the type of the atom to be created
         *                            i:e predicateNode, shemaNode
         * @return                 return -1 todo:// return void
         */
        Type atomese_argument(Handle&, const argument&, const id::procidure_types&);

        /**
         * Convert a combo argument to atomese.
         *
         * @param const builtin&       a builtin ref containing the combo builtin to be converted to atom
         * @param procidure_type&     set procidure type from the Link to be created
         *                            i:e predicate, shema
         * @return                 return Link type to be created
         */
        Type atomese_builtin(const builtin&, id::procidure_types&);
    }
}
#endif //MOSES_COMBO_ATOMESE_H
