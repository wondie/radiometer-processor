Requires Python 3.9 via Anaconda
GEOS 3.3 bin and dev or later. Install using OSGeo4W for Windows.

installation

conda install -c conda-forge pyinstaller

conda install -c conda-forge parallel

conda install -c conda-forge cmocean

conda install -c intel scikit-learn

conda install -c conda-forge astropy

conda install -c conda-forge georaster

conda install -c conda-forge cartopy

Inside  Anaconda Command Prompt  
pip install mgsub

pip install matplotlib==3.5.2

pip install pytexit

pip install sympy

pip install pytorch
# Go to anaconda shell:
cd .\libs\pytorch-minimize\
pip install -e .


pip install -r requirements.txt


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