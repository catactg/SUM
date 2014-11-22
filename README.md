
# Standardised unfold map (SUM) of the left atrium: regions definition for multimodal image analysis

The pipeline is split in parts.

run_std:              run standardisation on meshes generate by scard3d
                       (uses Meshlab and VMTK)
                       
run_currents:         currents registration
                       (uses deformetrica)
                       
run_sum:              make template disk (SUM)

run_quantification:   quantify circumferential coverage, longitudinal coverage and area coverage.


Dependencies:

python 2.7.6

vmtk http://www.vmtk.org

vtk 5.10.1 (included in vmtk)

numpy 1.8.0

Meshlab http://meshlab.sourceforge.net (skip if Poisson reconstruction is not desired)

deformetrica (http://www.deformetrica.org)


