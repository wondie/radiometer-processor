

Requires Python 3.9 via Anaconda, you can install the latest env and create an Python 3.9 env. 

conda create -n gdal_cf python=3.10
conda activate gdal_cf
conda config --env --add channels conda-forge

conda install gdal

GEOS 3.3 bin and dev or later. Install using OSGeo4W for Windows.

installation

conda install -c conda-forge pyinstaller

conda install -c conda-forge parallel

conda install -c conda-forge cmocean

conda install anaconda::scikit-learn

conda install -c conda-forge astropy

conda install -c conda-forge georaster

conda install -c conda-forge cartopy

conda install anaconda::pandas


conda install conda-forge::matplotlib

conda install anaconda::natsort

conda install anaconda::openpyxl

conda install conda-forge::xlsxwriter
conda install conda-forge::r-mgsub

conda install conda-forge::pytexit

conda install pytorch torchvision torchaudio pytorch-cuda=12.4 -c pytorch -c nvidia



To create portable exe file:
pyinstaller --clean .\RadiometerProcessorOne.spec

# Running the algorithm_generator.py
- To generate an algorithm Run `C:\Anaconda3\python.exe radiometer_processor/algorithm_generator.py`
- Choose the algorith with the best r-squared and errors and copy the output. 
- Get the output the equation in `radiometer_processor\data\result_[r-squared].png`

# Running the box_plot_line.py
- Run python `C:\Anaconda3\python.exe radiometer_processor/box_plot_line.py`
- Get the result `radiometer_processor\data\boxplot.png`

# Applying the SPM algorithm on UAS Imagery
- Run `\python.exe" outputProcessing/SPM_calculation.py`