#include <iostream>

#include <opencog/util/mt19937ar.h>
#include <opencog/util/Logger.h>
#include <moses/comboreduct/interpreter/eval.h>

#include <moses/moses/deme/deme_expander.h>
#include <moses/moses/metapopulation/metapopulation.h>

#include <moses/moses/representation/representation.h>
#include <moses/moses/moses/moses_main.h>
#include <moses/moses/optimization/optimization.h>
#include <moses/moses/scoring/scoring_base.h>

#include "pole_scoring.h"


using namespace std;
using namespace boost;
using namespace opencog;
using namespace moses;
using namespace reduct;


int main(int argc, char** argv)
{

    //set flag to print only cassert and other ERROR level logs on stdout
    opencog::logger().set_print_error_level_stdout();
    //read maximum evaluations and RNG seed from command line
    int max_evals;
    int seed;
    bool reduce=true;
    try {
        if(argc < 4){
            cerr << "Not enough arguments." << endl;
            throw std::length_error("Missing arguments.");
        }
        max_evals=lexical_cast<int>(argv[1]);
        seed=lexical_cast<int>(argv[2]);
        set_stepsize(1.25); //lexical_cast<double>(argv[3]));
        set_expansion(1.5); //lexical_cast<double>(argv[4]));
        set_depth(4) ; //exical_cast<int>(argv[5]));
        reduce = lexical_cast<bool>(argv[3]);
    } catch (...) {
        cerr << "Usage: " << argv[0] << " max_evals seed reduce?{0,1}" << endl <<
            "ann_combo_tree" << endl <<
            "Example:" << endl <<
            "- Arguments: 100 1 1" << endl <<
            "- Ann Combo Tree: ann($N1($I2 $I3 $I4 $I5 0.5 0.5 0.5 0.5) $N6($I2 $I3 $I4 $I5 0.5 0.5 0.5 0.5))" << endl <<
            "It uses 1.25 as step_size, 1.5 as expansion and 4 as depth." << endl;
        exit(1);
    }

    //read in seed tree
    combo_tree tr;
    cin >> tr;

    randGen().seed(seed);

    type_tree tt(id::lambda_type);
    tt.append_children(tt.begin(), id::ann_type, 1);

    //DOUBLE MARKOVIAN POLE TASK`
    const reduct::rule* si = &(ann_reduction());
    if(!reduce)
        si = &(clean_reduction());

    ann_pole2nv_bscore p2_bscore;
    behave_cscore cscorer(p2_bscore);

    univariate_optimization univ;
    deme_expander dex(tt, *si, *si, cscorer, univ);
    metapopulation metapop_pole2(tr, cscorer);

    boost::program_options::variables_map vm;
    jobs_t jobs;
    moses_parameters moses_param(vm, jobs, true, max_evals);
    moses_statistics st;
    run_moses(metapop_pole2, dex, moses_param, st);

    //change best combo tree back into ANN
    tree_transform trans;
    combo_tree best = metapop_pole2.best_tree();
    ann bestnet = trans.decodify_tree(best);

    //show best network
    cout << "Best network: " << endl;
    cout << &bestnet << endl;
    //write out in dot format
    bestnet.write_dot("best_nn.dot");

    CartPole *the_cart;
    the_cart = new CartPole(true,false);
    the_cart->nmarkov_long=true;
    the_cart->generalization_test=false;
    double fitness = the_cart->evalNet(&bestnet);
    delete the_cart;
    //for parameter sweep
    cout << fitness << endl;
}




