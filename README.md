# CoreMS LC Modules
A Python package for processing CoreMS assignments of ultrahigh mass resolution LC-ESI-MS data. 

### Functionality 
##### Functions that do not require CoreMS:
- Alignment of assigned features across a dataset
- Calculation of average assignment & feature parameters across a dataset:
    1. measured m/z
    2. calibrated m/z
    3. resolving power
    4. m/z error
    5. feature S/N
    6. confidence score 
- Gap filling of ambiguous assignments 
- Stoichiometric classifications 
- NOSC calucations 
- O/C, H/C, N/C calculations 
- Identification of significant assignment errors in a dataset, based on rolling average and standard deviation

##### Functions that require CoreMS:
- Determination of a feature's chromatographic dispersity
- Generation of calibrant list(s) for data calibration 
- QC checks of retention and intensity of an internal standard across a dataset 

### Installation

##### Local editable installation
    git clone https://github.com/deweycw/corems-tools.git corems-tools
    cd corems-tools/src
    python -m pip install -e . # from root directory of project  

##### Install pkg from TestPyPi
    python -m pip install --index-url https://test.pypi.org/simple/ --no-deps coremstools

### Build and upload to TestPyPi
See: https://packaging.python.org/en/latest/tutorials/packaging-projects/

##### Build prior to uploading to archives
    python -m pip install -U build
    python -m build

##### Upload to TestPyPi with Twine
    python -m pip install -U twine
    python -m twine upload --repository testpypi dist/*

### Example Workflow

##### Phase 1: Assignments & initial QC
1. Generate calibrants using calibrant generator function (requires CoreMS)
2. Perform assignments (requires CoreMS)
3. Calculate dispersity (requires CoreMS)
4. Evaluate retention and intensity of internal standard (requires CoreMS)
5. Generate assignment error plots (error v. m/z; error dist in each time window)

##### Phase 2: Alignment, gapfilling, blank correction, error flags
6. Align features across dataset
7. Gapfill across dataset 
8. Perform blank correction 
9. Flag features with potentially signficant assignment error (i.e., feature is flagged if difference between rolling average of assignment error (across dataset) and error of individual assignment exceeds 4x the standard deviation of the assignment error for the specific feature) 

##### Phase 3: Additional classifications 
10. Determine O/C, H/C, N/C & NOSC for features 
11. Determine stoichiometric classifcations



