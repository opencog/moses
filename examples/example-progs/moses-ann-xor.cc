
#include <iostream>

#include <opencog/util/mt19937ar.h>
#include <opencog/util/Logger.h>

#include <moses/comboreduct/interpreter/eval.h>

#include <moses/moses/deme/deme_expander.h>
#include <moses/moses/metapopulation/metapopulation.h>

#include <moses/moses/moses/moses_main.h>
#include <moses/moses/scoring/scoring_base.h>
#include <moses/moses/optimization/optimization.h>
#include <moses/moses/representation/representation.h>

#include "ann_xor_scoring.h"

using namespace moses;
using namespace reduct;
using namespace boost;
using namespace std;
using namespace opencog;


int main(int argc, char** argv)
{

    //set flag to print only cassert and other ERROR level logs on stdout
    opencog::logger().set_print_error_level_stdout();

    //read in maximum evaluations and RNG seed from command line
    int max_evals;
    int seed;
    bool reduce=true;
    try {
        if(argc < 7){
            cerr << "Not enough arguments." << endl;
            throw std::length_error("Missing arguments.");
        }
        max_evals=lexical_cast<int>(argv[1]);
        seed=lexical_cast<int>(argv[2]);
        set_stepsize(lexical_cast<double>(argv[3]));
        set_expansion(lexical_cast<double>(argv[4]));
        set_depth(lexical_cast<int>(argv[5]));
        reduce = lexical_cast<bool>(argv[6]);
    } catch (...) {
        cerr << "Usage: " << argv[0] << "max_eval seed step_size expansion depth reduce?{0,1}" << endl <<
            "ann_combo_tree" << endl <<
            "Example:" << endl <<
            "- Arguments: 100 1 1 2 5 1" << endl <<
            "- Ann Combo Tree: ann($N1($I2 $I3 $I4 0.0 0.0 0.0))" << endl;
        exit(1);
    }

    // read in the seed combo tree from stdin
    // a default seed is provided in file xor_tree
    combo_tree tr;
    cin >> tr;

    randGen().seed(seed);

    // this will let representation building know that we are dealing
    // with an ANN, what it says is that the type of combo_tree to
    // evolve are ANN generators (no input, output an ANN)
    type_tree tt(id::lambda_type);
    tt.append_children(tt.begin(), id::ann_type, 1);


    // binary XOR task
    ann_xor_bscore bscore;
    behave_cscore cscorer(bscore);

    const reduct::rule* si = &(ann_reduction());
    if(!reduce)
        si = &(clean_reduction());

    univariate_optimization univ;
    deme_expander dex(tt, *si, *si, cscorer, univ);
    metapopulation metapop(tr, cscorer);

    boost::program_options::variables_map vm;
    jobs_t jobs;
    moses_parameters moses_param(vm, jobs, true, max_evals);
    moses_statistics st;
    run_moses(metapop, dex, moses_param, st);

    //transform the best combo tree into an ANN
    tree_transform trans;
    combo_tree best = metapop.best_tree();
    ann bestnet = trans.decodify_tree(best);


    //look at xor outputs
    double inputs[4][3] = { {0.0, 0.0, 1.0},
                            {0.0, 1.0, 1.0},
                            {1.0, 0.0, 1.0},
                            {1.0, 1.0, 1.0}};

    int depth = bestnet.feedforward_depth();

    for (int pattern = 0;pattern < 4;++pattern) {
        bestnet.load_inputs(inputs[pattern]);
        for (int x = 0;x < depth;++x)
            bestnet.propagate();
        cout << "Input [ " << inputs[pattern][0] << " " << inputs[pattern][1] << " ] : Output " << bestnet.outputs[0]->activation << endl;

    }

    //save the best network (can be viewed in any dot viewer)
    cout << "Best network: " << endl;
    cout << &bestnet << endl;
    bestnet.write_dot("best_nn.dot");

    //write out the best score (to be used in parameter sweep)
    cout << metapop.best_score() << endl;
}
