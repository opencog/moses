ADD_SUBDIRECTORY (moses)

# direct file and path creation for stuff that doesn't need to be compiled
file(MAKE_DIRECTORY moses)
# module init
file(COPY moses/__init__.py DESTINATION moses)
##
