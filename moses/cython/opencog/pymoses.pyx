__author__ = 'Cosmo Harrigan'

from cpython.version cimport PY_MAJOR_VERSION
from libc.stdlib cimport malloc, free
import shlex
import tempfile
import csv
import sys

class MosesException(Exception):
    pass

class MosesCandidate(object):
    def __init__(self, score = None, program = None, program_type = None):
        self.score = score
        self.program = program
        self.program_type = program_type

    def eval(self, arglist):
        if self.program_type != "python":
            raise MosesException('Error: eval method is only defined for '
                                 'candidates with program_type of python.')
        if len(arglist) == 0:
            raise MosesException('Error: eval method requires a list of input '
                                 'values.')

        namespace = {}
        exec(self.program, namespace)
        return namespace.get('moses_eval')(arglist)

cdef class moses:
    def run(self, input = None, args = "", python = False , scheme = False):
        """
        Invokes MOSES in supervised learning mode to learn candidate solutions
        for a given training set.

        Parameters:
            input (list of lists) - training data for regression [optional]
                Example: input=[[0, 0, 0], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
            args (string) - arguments for MOSES (see MOSES documentation)
                            [optional]
            python (bool) - if True, return Python instead of Combo [optional]
        Either input or args must be provided, otherwise, MOSES would have no
        input

        Output:
        Returns a collection of candidates as MosesCandidate objects.
        Each MosesCandidate contains:
            score (int)
            program (string)
            program_type (string - Enumeration: python, combo)
            eval (Runnable method - Only valid for program_type python) Run
            this method to evaluate the model on new data.
        """
        if (input is None or input == "") and args == "":
            raise MosesException('Error: input and args cannot both be empty. '
                                 'You should pass your input using the input '
                                 'parameter, or refer to it in the args '
                                 'parameter.')

        # Create temporary files for sending input/output to the moses_exec
        # function
        if input is not None:
            if PY_MAJOR_VERSION < 3:
                input_file = tempfile.NamedTemporaryFile()
            else:
                input_file = tempfile.NamedTemporaryFile(mode='w+', encoding='utf-8')
            input_file_builder = csv.writer(input_file, delimiter = ',')
            input_file_builder.writerows(input)
            input_file.flush()

        output_file = tempfile.NamedTemporaryFile()

        # Process the argument list for moses_exec
        _args_list = []

        if input is not None:
            _args_list.extend(['-i', input_file.name])
        if python and not scheme :
            _args_list.extend(['--output-format', 'python'])
        if scheme and not python :
            _args_list.extend(['--output-format', 'scheme'])


        _args_list.extend(['-o', output_file.name])
        _args_list.extend(shlex.split(args))

        self._run_args_list(_args_list)
        # Process the output file
        output = output_file.file.read()

        candidates = []
        # Python header declared in moses/moses/moses/types.h
        # (ostream_combo_tree_composite_pbscore_python)
        python_header = b"#!/usr/bin/env python"

        if len(output) == 0:
            raise MosesException('Error: No output file was obtained from '
                                 'MOSES. Check to make sure the input file '
                                 'and arguments provided were valid.')

        # Python output
        elif output.splitlines()[0].startswith(python_header):
            output_list = [python_header + b"\n" + element for
                           element in output.split(python_header)[1:]]

            for candidate in output_list:
                program = candidate
                if b"#score: " in program:
                    score = int(program.split(b"#score: ")[1].splitlines()[0])
                else:
                    raise MosesException('Error: A score value was expected '
                                         'but not found in the Python '
                                         'program.')
                candidates.append(MosesCandidate(score = score,
                                                 program = program,
                                                 program_type = "python"))

        

        

       # Combo or scheme output
        else:
            output_list = [element for element in output.splitlines()]

            for candidate in output_list:
                # XXX FIXME, this is an utterly horrible way to get
                # MOSES output, because it tries to parse what MOSES
                # prints, and that will change for abitrary reasons.
                # The output string is NOT a part of the formal API
                # for moses!
                score = int(candidate.partition(b' ')[0])
                rest = candidate.partition(b' ')[2]
                # weight = int(rest.partition(b' ')[0])
                # rest = rest.partition(b' ')[2]
                program = rest.partition(b'[')[0]
                if scheme :

                    candidates.append(MosesCandidate(score = score,
                                                    program = program,
                                                    program_type = "scheme"))
                else :

                    candidates.append(MosesCandidate(score = score,
                                                    program = program,
                                                    program_type = "combo"))

        return candidates
   
    
    def write_scheme(self, candidates=[]):
        '''
        writes out programs with the best scores 
        on a separate scheme file so it can be used for other purposes 
        (ie - it can be loaded into opencog atomspace)

        '''
        
        if len(candidates) == 0:
            raise MosesException('Error: write_scheme method requires a list of input '
                                 'values.')
        if candidates[0].program_type != "scheme":
            raise MosesException('Error: eval method is not defined for '
                                 'candidates with program_type other than scheme.')
        else:
            best_score = candidates[0].score
            for candidate in candidates:
                if candidate.score > best_score:
                    best_score = candidate.score
            output_file = "moses_result.scm"
            of = open(output_file,"w")
            for candidate in candidates:
                if candidate.score == best_score:
                    of.write(candidate.program + "\n")




    def run_manually(self, args=""):
        """
        Invokes MOSES without any extra integration.
        Useful for interacting with MOSES non-programatically via stdout.
        """
        if args == "":
            raise MosesException('Error: No arguments provided')
        self._run_args_list(shlex.split(args))

    def _run_args_list(self, args_list):
        args_list.insert(0, "moses")
        cdef char **c_argv
        if PY_MAJOR_VERSION < 3:
            args_list = [bytes(x.decode("utf8")) for x in args_list]
        else:
            args_list = [bytes(x, "utf8") for x in args_list]
        c_argv = <char**>malloc(sizeof(char*) * len(args_list))
        for idx, s in enumerate(args_list):
            c_argv[idx] = s
        try:
            moses_exec(len(args_list), c_argv)
        except RuntimeError, ex:
            raise MosesException('Error: exception occurred calling C++ MOSES.')

            # WTF. FIXME. Something about exceptions is borken.
            # AttributeError: 'RuntimeError' object has no attribute 'message'

            # raise MosesException('Error: exception occurred when calling C++ '
            #                     'MOSES. Exception message:\n' + ex.message)
        finally:
            free(c_argv)
