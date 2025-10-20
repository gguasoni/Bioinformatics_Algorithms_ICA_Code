## Task
The focus of this assignment was two-fold:
1. Build a command line interface (CLI) for the Blast101 programme. The CLI had to...
   - intercept (incorrect) inputs early.
   - ask the user to enter query sequence and database, whilst keeping all other parameters constant.
   - check that the sequence and database are valid
   - print the state of the input variables before running the code
   - allow the user the opportunity to swap between unning a BLAST search, running a Smith-Waterman search, or running the code for the statistical analysis.
   - ask the user if it wants to run the tests against the BLAST algorithm and CLI.

2. Writing tests for the Blast101 application.

## Script Details
Students were given a folder of scripts implementing BLAST101 and Smith-Waterman algorithms to work around and through to establish the CLI and write the application tests. As I am not the author of this code, I am not including it in this repository.

### Authored scripts
- The file to execute the script is final_ba_script.py
- The file containing extra functions for the command_line interface is ui_script.py
- test_script.py contains the functions for testing the command line interface, the BLAST analysis, and the Smith-Waterman analysis. 

## How I would improve this assignment

## Coding References
The following websites were consulted at various points during the completion of this assignment, and the code they contain were necessarily modified for the assignment. 
https://www.geeksforgeeks.org/command-line-interface-programming-python/
https://docs.python.org/3/library/argparse.html#the-add-argument-method
https://docs.python.org/3/tutorial/datastructures.html#sets
https://tiefenauer.github.io/blog/smith-waterman/
https://docs.python.org/3/reference/datamodel.html#function.__name__
https://docs.python.org/3/tutorial/errors.html

