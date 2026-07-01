# CMake Build Notes

The project uses CMake presets for configuration and build.

## Requirements

You need:

- `cmake`
- `ninja`
- one C compiler:
  - MSVC from Visual Studio Build Tools
  - or GCC from MinGW/MSYS2
  - or Clang/LLVM

## Presets

- `default`: Release build with `Ninja`
- `debug`: Debug build with `Ninja`
- `mingw`: Release build with `MinGW Makefiles`

## Build

The project path contains Cyrillic characters. MinGW cannot compile sources
from such a path directly, so use the wrapper that creates a temporary ASCII
junction and adds the MinGW DLL directory to `PATH`:

```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1 -Fresh
```

Subsequent builds:

```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1
```

## Windows Examples

### MSVC

Run the commands from `x64 Native Tools Command Prompt for VS` or from a PowerShell session where the MSVC environment is loaded.

```powershell
cmake --preset default
cmake --build --preset default
```

### MSYS2 / MinGW (manual alternative)

```powershell
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake --preset mingw
cmake --build --preset mingw
```

## Report Checklist

The repository contains program sources only; prepare the lab report as a
separate document according to the TUSUR standard and the UMP structure.

Recommended report structure:

- title page;
- task statement and variant data;
- brief theory for each method;
- algorithm descriptions or block diagrams;
- program listing or key fragments;
- calculation results with intermediate iterations;
- accuracy and residual estimates in scientific notation;
- conclusion.
