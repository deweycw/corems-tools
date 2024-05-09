# CoreMS LC Modules
A Python package for processing CoreMS assignments of ultrahigh mass resolution LC-ESI-MS data. 

#### Functionality 
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
- NOSC calculations 
- O/C, H/C, N/C calculations 
- Identification of significant assignment errors in a dataset, based on rolling average and standard deviation
- Identification of features in samples that also appear in blanks 


##### Functions that require CoreMS:
- Determination of a feature's chromatographic dispersity
- Generation of calibrant list(s) for data calibration 
- QC checks of retention and intensity of an internal standard across a dataset 

#### Build/install commands:
See: https://packaging.python.org/en/latest/tutorials/packaging-projects/

###### For editable local installation 
    python -m pip install -e . # from root directory of project  
        
###### Build prior to uploading to archives
    python -m pip install -U build
    python -m build

###### Upload to TestPyPi with Twine
    python -m pip install -U twine
    python -m twine upload --repository testpypi dist/*

###### Install pkg from TestPyPi
    python -m pip install --index-url https://test.pypi.org/simple/ --no-deps PACKAGE_NAME



