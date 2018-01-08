# Derivative-free Global Optimization Using Space-filling Curves
## DIRECT
Python implementation of Dividing Rectangles global search algorithm
- based on: http://www4.ncsu.edu/~ctk/Finkel_Direct
- [insert algo/psuedo code]

## Hilbert Curve
C++ implementation of Space-filling curve ...
- based on: [John Skilling] [link]
- [insert algo/psuedo code]


## Requirements
- Python(CPython)3.4+
- numpy
- appropriate C++ compiler (VS,g++,gcc, etc.)

## Setup
MS Visual Studio (2013 and above) settings:
1. Go to `Project > Hilbert Properties... > Configuration Properties`. With All Configurations:
   - Chose VC++ Directories setting. Append `$(PYTHONPATH)\include` to Include Directories for the Python header file `<Python.h>` and append `$(PYTHONPATH)\libs` to Library Directories for the linker library file `pythonXX.lib` or `pythonXX_d.lib`.
   - Choose General settings and in Project Defaults, change Configuration Type to 'Dynamic Library (.dll)'.
   - Under General settings, change Target Extension to '.pyd'. 
   <br /><br />
   Only for Debug Configuration:
   - Under General settings, change Target Name to `$(ProjectName)_d`.

2. Chose either Release or Debug under Solution Configurations. Note that for Debug mode, you need to have Python debug binaries downloaded first (requires VS 2015 or later).

3. With Release configuration, Build Visual Studio project with target build directory as `DIRECT\build`. Run `dir DIRECT\build` to locate `Hilbert.pyd` file.
<br />
For Debug mode, steps are similar as above. `Hilbert_d.pyd` file will be created under `Hilbert\Debug` and you may use the `pythonXX_d.exe` interpreter instead.

## Run
[how to call/use wrapper]
CMD:
```Shell
python DIRECT\src\main.py	# invoke Direct.run()
python DIRECT\test.py		# invoke Hilbert.lib
```

## Project Structure
Main project files:
```
root
|
|- src
|	|
|	|- direct.py
|	|- helper.py
|	|- main.py
|
|- build
|	|
|	|- Hilbert.lib
|	|- Hilbert.pyd
|
|- test.py
```
[file contents gist]

NOTE: README in progress.
