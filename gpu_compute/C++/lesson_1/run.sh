326  g++ -o hello main.cpp
  327  ./hello

ctrl shift b
./main

# Line 1: #include <iostream>
# Purpose: Preprocessor directive to include the Input/Output Stream library
# What it does: Makes C++ input/output functionality available (like cout, cin)
# Analogy: Like importing a module in Python or including a header in C

# Line 2: using namespace std;
# Purpose: Brings the standard namespace into scope
# What it does: Allows you to use cout instead of std::cout, endl instead of std::endl
# Note: While convenient for learning, in larger projects it's better to use std::cout explicitly to avoid naming conflicts

# Line 4: int main() {
# Purpose: Defines the main function - the entry point of every C++ program
# Return type: int (integer) - the program returns an exit status to the operating system
# Parameters: Empty parentheses () mean no command-line arguments are expected

# Line 5: cout << "Hello, world!" << endl;
# cout: Standard output stream object (console output)
# <<: Stream insertion operator - "sends" data to the output stream
# "Hello, world!": String literal to be displayed
# endl: End line - inserts a newline and flushes the output buffer
# Equivalent to: cout << "Hello, world!\n"; (but endl also flushes the buffer)

# Line 6: return 0;
# Purpose: Returns exit status to the operating system
# 0: Conventionally means "successful execution"
# Non-zero values: Typically indicate error conditions

# Line 7: }
# Purpose: Closes the main function block
# Key C++ Concepts Demonstrated:
# Header includes (#include)
# Namespaces (using namespace)
# Function declaration (int main())
# Stream I/O (cout, <<)
# Return statements and program exit codes
# This is the canonical "Hello, World!" program that introduces fundamental C++ syntax and structure!