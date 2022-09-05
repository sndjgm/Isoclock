# Isoclock
A novel off-line software for data reduction of LA-ICP-MS U-Pb dating with reference containing variable common Pb
Isoclock is designed to deduct common Pb for common-Pb-bearing materials and calibrate the U-Pb isotope fractionation. The software contains several processing steps for the raw files of mass spectrometry, such as data import and view, background correction and filtering of outliers, calculation of common Pb for reference materials, fractionation calibration, and age calculation. Isoclock is written using the free and open-source Python language and can either run the code directly or operate using a graphical interface. It is compatible with raw data files from widely used modern ICP-MS instruments and allows for extended data interfaces.



I. Software operation procedures

1.1 Start-up and functional block division

The following steps 1.1.1-1.1.5 is run before Isoclock runs for the first time.

1.1.1	Python 3.9 is necessary to run the code. Download from https://www.python.org/downloads/ and follow the installation.

1.1.2 Download or clone this repository.

1.1.3 Open terminal/cmd and navigate to the Isoclock folder.

 cd path/to/folder/ Isoclock

1.1.4 Instal python libraries required for Isoclock.

	pip install -r requirements.txt

1.1.5 Run Isoclock from python.

	python Isoclock.py

If everything is already installed, follow only steps1.1.5. If you are Windows user, you can also run the Isoclock.exe directly.



