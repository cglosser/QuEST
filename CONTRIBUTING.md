# How to contribute

First of all, thank you for taking time to check out QuEST (and perhaps the associated articles). Since you're here it means you have at least a passing interest in contributing some work to QuEST and we'd really love to have both your expertise and help. Science only moves forward through collaboration. 

First thing's first, if you haven't already, please send the project maintainer (@cglosser) an email or open an issue introducing yourself and describing your (ideas for) changes. Until the project gets a little more momentum that will really help mitigate conflicting chages.

## Testing

Implementing tests is perhaps the easiest and most important way you can help develop QuEST. As this tool is fundamentally mathematical in nature, there are zillions of routines to test and benchmark against fully analytical (pencil & paper) solutions. QuEST leverages both `boost.test` and [`rapidcheck`](https://github.com/emil-e/rapidcheck) for testing; please use the former and consider using the latter for testing where appropriate. If you do work on something mathematical, please also be sure to include a copy of the necessary formulae and/or expressions for what you've done.

## Submitting changes

To submit your changes, please send a GitHub Pull Request _to the `development` branch of QuEST_ with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)). We can always use more test coverage. Please follow our coding conventions (below) and make sure all of your commits are atomic (one feature per commit).

Always write a clear log message for your commits. One-line messages are fine for small changes, but bigger changes should look like this:

    $ git commit -m "A brief summary of the commit
    > 
    > A paragraph describing what changed and its impact."

We will love you forever if you refer to the issue numbers here on GitHub that relate to your changes (e.g. if you fix a bug described in issue #37, mention #37 somewhere in your commit).

## Coding conventions

Start reading our code and you'll get the hang of it. For the most part, we optimize for readability. If you want to submit changes that bring some code more in-line with these conventions, please do so!

  * The top-level `.clang-format` file contains the preferred styling rules for all C++ code. Please follow that and run that on your source files judiciously. 
  * ALWAYS put spaces after list items and method parameters (`[1, 2, 3]`, not `[1,2,3]`), around operators (`x += 1`, not `x+=1`), and around hash arrows.
  * Set class names in `CamelCase` and variable/function names in `snake_case`.
  * Short variable names are right out. Be descriptive and try to make your code read like a sentence. `chi[i][j]` doesn't work nearly as well as `damping_parameters[particle_idx][derivative_idx]`.
  * This is open source software. Consider the people who will read your code, and make it look nice for them. It's sort of like driving a car: Perhaps you love doing donuts when you're alone, but with passengers the goal is to make the ride as smooth as possible.
  
### References and Pointers

   * Prefer passing references/values instead of pointers unless a given function _needs_ a pointer (likely to give ownership of a shared resource or abstract type). This pushes the responsibility of checking for null pointers higher in the function tree, thereby reducing error-checking boilerplate.

## Academic integrity

QuEST is an open-source platform for research and a lot of development work went into it. Please make use of the tool as you need, but also please consider paying it back with even small features or fixes---every little bit greatly helps the project! And if you do use QuEST in your own research, please don't heasitate to let us know; we'd be thrilled to put a link to your work on the main page!

Thanks,
Connor Glosser, Active QuEST Maintainer (guildmaster?)
