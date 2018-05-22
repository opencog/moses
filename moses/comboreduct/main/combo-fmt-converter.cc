/*
 * moses/comboreduct/main/combo-fmt-converter.cc
 *
 * Copyright (C) 2015 OpenCog Foundation
 * All Rights Reserved
 *
 * Written by Nil Geisweiller
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

#include <fstream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <opencog/atoms/base/Handle.h>
#include <opencog/atoms/base/Atom.h>
#include <moses/comboreduct/combo/combo_atomese.h>
#include <opencog/atoms/base/Link.h>
#include <opencog/atoms/base/Node.h>

#include "../combo/iostream_combo.h"
#include "../type_checker/type_tree.h"

using namespace boost::program_options;
using namespace boost::algorithm;
using namespace std;
using namespace opencog;
using namespace combo;
using boost::trim;

// using namespace ant_combo;

// Structure containing the options for the combo-fmt-converter program
struct pgrParameters
{
	vector<string> combo_programs;
	vector<string> combo_programs_files;
	string log_level;
	string output_file;
	string output_format_str;
};

/**
 * Read and parse the combo-fmt-converter program arguments.
 * Return the parsed results in the parameters struct.
 */
pgrParameters parse_program_args(int argc, char** argv)
{
	// program options, see options_description below for their meaning
	pgrParameters pa;

	// Declare the supported options.
	options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Produce help message.\n")

		("combo-program,c",
		 value<vector<string>>(&pa.combo_programs),
		 "Alternative way to entering the programs to convert "
		 "(as opposed to via stdin, not that using this option will "
		 "disable reading on the stdin). "
		 "Remember to escape the $ signs or to wrap them in single "
		 "quotes as to not let the shell perform variable expansions "
		 "on the combo program variables."
		 "This option may be "
		 "used several times to convert several programs.\n")

		("log-level,l",
         value<string>(&pa.log_level)->default_value("INFO"),
         "Log level, possible levels are NONE, ERROR, WARN, INFO, "
         "DEBUG, FINE. Case does not matter.\n")

		("combo-programs-file,C",
		 value<vector<string>>(&pa.combo_programs_files),
		 "Other alternative way (disables stdin). "
		 "Indicate the path of a file containing combo "
		 "programs (seperated by newlines) to convert. "
		 "This option may be used several times to progress several "
		 "files.\n")

		("output-file,o",
		 value<string>(&pa.output_file),
		 "File where to save the results. If empty then outputs on "
		 "stdout.\n")

		("output-format,f",
		 value<string>(&pa.output_format_str)->default_value("combo"),
		 "Supported output formats are combo, python and scheme.\n")

		;

	variables_map vm;
	store(parse_command_line(argc, argv, desc), vm);
	notify(vm);

	if (vm.count("help") || argc == 1) {
		cout << desc << "\n";
		exit(1);
	}

	// Remove old log_file before setting the new one.
	const string log_filename = "combo-fmt-converter.log";
	remove(log_filename.c_str());
	logger().set_filename(log_filename);
	trim(pa.log_level);
	Logger::Level level = logger().get_level_from_string(pa.log_level);
	if (level != Logger::BAD_LEVEL)
		logger().set_level(level);
	else {
		cerr << "Error: Log level " << pa.log_level
		     << " is incorrect (see --help)." << endl;
		exit(1);
	}
	logger().set_backtrace_level(Logger::ERROR);
	logger().set_print_error_level_stdout();

	return pa;
}

vector<string> get_all_combo_tree_str(const pgrParameters& pa)
{
    vector<string> res(pa.combo_programs);     // from command line

    // from files
    for (const string& combo_programs_file : pa.combo_programs_files) {
        ifstream in(combo_programs_file.c_str());
        if (in) {
            while (in.good()) {
                string line;
                getline(in, line);
                if (line.empty())
                    continue;
                res.push_back(line);
            }
        } else {
            logger().error("Error: file %s can not be found.",
                           combo_programs_file.c_str());
            exit(1);
        }
    }

    return res;
}
// Convert a Combo_tree
Handle convert(const combo_tree &tr){
    return atomese_combo(tr);
}
// Convert a single combo program gotten from istream
void convert(istream& in, ostream& out,
             const pgrParameters& pa,
             combo::output_format fmt) {
	// Parse combo tree
	string line;
	getline(in, line);
	if (line.empty())
		return;

	// Remove useless wrapping double quotes
	if (starts_with(line, "\"") and ends_with(line, "\""))
		line = line.substr(1, line.size() - 2);

	vector<string> variables = parse_combo_variables(line);
	combo_tree tr = str2combo_tree(line, variables);
	type_tree tt = infer_type_tree(tr);

	// Write to the right format
	ostream_combo_tree(out, tr, variables, fmt) << std::endl;
}

int main(int argc, char** argv)
{
	// Parse program options
	pgrParameters pa = parse_program_args(argc, argv);
	
	// Parse output format
	combo::output_format fmt = parse_output_format(pa.output_format_str);

    // read combo program strings given by an alternative way than
	// stdin
	vector<string> combo_trees_str = get_all_combo_tree_str(pa);

	ostream* out = pa.output_file.empty()?
		&cout :	new ofstream(pa.output_file.c_str());

	if (combo_trees_str.empty()) {
		// If empty then use stdin
		while (cin.good())
			convert(cin, *out, pa, fmt);
	} else {
		// Don't use stdin
		for (const string& str : combo_trees_str) {
			stringstream ss(str);
			convert(ss, *out, pa, fmt);
		}
	}

	if (not pa.output_file.empty())
		delete out;

	return 0;
}
