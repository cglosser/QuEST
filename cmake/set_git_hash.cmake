# Adapted from Irshad Pananilath -- http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/

function(set_git_macros target)
  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    FIND_PACKAGE(Git)
    IF(GIT_FOUND)
      # Get the current working branch
      execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

      # Get the latest abbreviated commit hash of the working branch
      execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

      target_compile_definitions(${target} PUBLIC 
        __GIT_HASH__="${GIT_COMMIT_HASH}"
        __GIT_BRANCH__="${GIT_BRANCH}"
      )
    ELSE(GIT_FOUND)
      target_compile_definitions(${target} PUBLIC 
        __GIT_HASH__="UNAVAILABLE"
        __GIT_BRANCH__="UNAVAILABLE"
      )
    ENDIF(GIT_FOUND)
  ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
endfunction()
