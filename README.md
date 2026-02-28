# My C Project

This project is a simple C application that demonstrates the basic structure of a C program. It includes a main function that prints "Hello, World!" to the console.

## Project Structure

```
my-c-project
├── src
│   └── main.c        # Contains the main function
├── include
│   └── main.h        # Header file for function declarations
├── Makefile           # Build instructions
├── .vscode
│   └── tasks.json     # VS Code tasks configuration
├── .gitignore         # Git ignore file
└── README.md          # Project documentation
```

## Building the Project

To build the project, you can use the provided Makefile. Open a terminal in the project directory and run:

```
make
```

This will compile the source files and create the executable.

## Running the Application

After building the project, you can run the application with the following command:

```
./main
```

This will execute the program and print "Hello, World!" to the console.

## Contributing

Feel free to contribute to this project by submitting issues or pull requests.

## Build
```PowerShell
rm -r build
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake --preset default
cmake --build --preset default
```