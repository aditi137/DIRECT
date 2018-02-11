# Derivative-free Global Optimization Using Space-filling Curves
## DIRECT
Python implementation of Dividing Rectangles global search algorithm
- based on: DIRECT Optimization Algorithm User Guide, Dan Finkel <br />
  http://www4.ncsu.edu/~ctk/Finkel_Direct/DirectUserGuide_pdf.pdf

## Hilbert Curve
C++ implementation of Space-filling curve ...
- based on: Programming the Hilbert curve, John Skilling <br />
  http://ratml.org/misc/hilbert_curve.pdf

## Requirements
- Python(CPython)3.4+
- numpy
- pytest

## Setup 
Python Project IDE Settings:
1. Add `DIRECT\src` to the Project Source Path.

2. Set the Python Interpreter path to your `python.exe` installation directory and the Test Runner for the Project as Py.test.

## Run
```Shell
python DIRECT\src\main.py	# invoke Direct.run()
```

## Project Structure
Main project files:
```
root
|
|- src
|	|
|	|- _hilbert.py
|	|- direct.py
|	|- helper.py
|	|- main.py
|
|- test_Hilbert.py
```
[file contents gist]

NOTE: README in progress.
