ADD_SUBDIRECTORY(combo)
ADD_SUBDIRECTORY(interpreter)
ADD_SUBDIRECTORY(type_checker)
ADD_SUBDIRECTORY(table)
ADD_SUBDIRECTORY(reduct)
ADD_SUBDIRECTORY(ant_combo_vocabulary)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(man)
ADD_SUBDIRECTORY(crutil)

ADD_LIBRARY(comboreduct SHARED
	crutil/exception.cc

	combo/action.cc
	combo/action_symbol.cc
	combo/argument.h
	combo/assumption.cc
	combo/builtin_action.cc
	combo/convert_ann_combo.h
	combo/definite_object.cc
	combo/descriptions.cc
	combo/enum_type.cc
	combo/indefinite_object.cc
	combo/iostream_combo.cc
	combo/message.cc
	combo/perception.cc
	combo/procedure_call.cc
	combo/procedure_repository.cc
	combo/similarity.cc
	combo/simple_nn.cc
	combo/vertex.cc

	table/table.cc
	table/table_io.cc

	interpreter/interpreter.cc
	interpreter/eval.cc

	type_checker/type_tree_def.cc
	type_checker/type_tree.cc

	reduct/action_reduction.cc
	reduct/action_rules.cc
	reduct/branch_rules.cc
	reduct/clean_reduction.cc
	reduct/contin_reduction.cc
	reduct/contin_rules.cc
	reduct/flat_normal_form.cc
	reduct/full_reduction.cc
	reduct/general_rules.cc
	reduct/logical_reduction.cc
	reduct/ann_reduction.cc
	reduct/logical_rules.cc
	reduct/fold_reduction.cc
	reduct/fold_rules.cc
	reduct/meta_rules.cc
	reduct/mixed_reduction.cc
	reduct/mixed_rules.cc
	reduct/perception_reduction.cc
	reduct/perception_rules.cc
	reduct/reduct.cc
)

TARGET_LINK_LIBRARIES (comboreduct
	${COGUTIL_LIBRARY}
	${Boost_DATE_TIME_LIBRARY}
	${Boost_THREAD_LIBRARY}
	${Boost_REGEX_LIBRARY}
)

ADD_LIBRARY(comboant SHARED
	ant_combo_vocabulary/ant_builtin_action.cc
	ant_combo_vocabulary/ant_perception.cc
	ant_combo_vocabulary/ant_action_symbol.cc
	ant_combo_vocabulary/ant_indefinite_object.cc
	ant_combo_vocabulary/ant_combo_vocabulary.cc
)

TARGET_LINK_LIBRARIES (comboant
	comboreduct
	${COGUTIL_LIBRARY}
)


#install library
IF (WIN32)
	INSTALL(TARGETS comboreduct comboant
                EXPORT MosesTargets
                DESTINATION "lib${LIB_DIR_SUFFIX}/moses")
ELSE (WIN32)
	INSTALL(TARGETS comboreduct comboant
                EXPORT MosesTargets
           LIBRARY DESTINATION "lib${LIB_DIR_SUFFIX}"   # lib*.so files
           ARCHIVE DESTINATION "lib${LIB_DIR_SUFFIX}")  # lib*.a files
ENDIF (WIN32)
