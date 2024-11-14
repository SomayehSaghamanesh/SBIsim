# SBIsim

A software to simulate X-ray phase-contrast and dark-field imaging techniques 
(Speckle-Based Imaging (SBI) and Propagation-Based imaging (PBI) for now). 
Below options are available in version 1.0: 
- monochromatic or polychromatic X-ray sources
- parallel (not implemented yet) or cone-beam X-ray beams/geometries
- Thin/thick samples
- different object and/or diffuser materials
** virtual object has not yet been implemented in version 1.0.

# How to Run
The software doesn't need any installation. Download only the SBIsim-1.0 file (only for Linux). 
Then follow below steps:
1- In linux, open the terminal (ctrl+shift+T), then:
2- $ cd /path/to/your_downloaded_file
3- $ sudo chmod +x SBIsim-1.0
4- $ ./SBIsim-1.0
above will run the software as well as show some verboses inside the terminal while running the software.

** As an alternative, you can just skip step 4 and click on the application file (SBIsim-1.0), if you already have 
FUSE and AppImageLauncher installed on your linux. Those utilities basically make SBIsim-1.0 as an
usual application installed on (and integrated in) your linux system, where you can call it from the applications
list and run it by one click. However, no verbose will be shown by this way and the user will not be informed
when the simulation is finished. In next versions, a status bar will be added for the processing.  
