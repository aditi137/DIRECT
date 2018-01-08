# Derivative-free Global Optimization Using Space-filling Curves
## DIRECT
Python implementation of Dividing Rectangles global search algorithm
- based on: http://www4.ncsu.edu/~ctk/Finkel_Direct
- [insert algo/psuedo code]

## Requirements
- Python(CPython)3.4+
- numpy

## Setup and Build
CMD:
Execute the following command based on the version of Visual Studio installed:
    Visual Studio 2010 (VS10): `SET VS90COMNTOOLS=%VS100COMNTOOLS%`
    Visual Studio 2012 (VS11): `SET VS90COMNTOOLS=%VS110COMNTOOLS%`
    Visual Studio 2013 (VS12): `SET VS90COMNTOOLS=%VS120COMNTOOLS%`
    Visual Studio 2015 (VS14): `SET VS90COMNTOOLS=%VS140COMNTOOLS%`

```Shell
SET DISTUTILS_USE_SDK=1
SET MSSdk=1
python setup.py build
```

## Run
[how to call/use wrapper]
CMD:
```Shell
python src\main\main.py
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
|- setup.py
```
[file contents gist]

NOTE: README in progress.