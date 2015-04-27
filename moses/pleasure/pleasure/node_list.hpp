#include <set>
#include <istream>

#include <moses/comboreduct/combo/vertex.h>
#include <moses/comboreduct/ant_combo_vocabulary/ant_combo_vocabulary.h>
#include <moses/comboreduct/combo/builtin_action.h>
#include <moses/comboreduct/combo/perception.h>
#include <moses/comboreduct/combo/action_symbol.h>
#include <moses/comboreduct/combo/indefinite_object.h>

namespace pleasure {
    typedef std::set<moses::combo::vertex> node_list;

    //std::istream& stream_to_node_list(std::istream& is, node_list& list);
    void stream_to_node_list(std::istream& is, node_list& list);
}

