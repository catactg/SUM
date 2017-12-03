# SUM: Standardized unfold map


Author: Catalina Tobon Gomez (catactg@gmail.com)

## About


This code generates a *standardized unfold map* as described in

> Standardized unfold mapping: a technique to permit left atrial regional data display and analysis.
> Williams SE, Tobon-Gomez C, Zuluaga MA, Chubb H, Butakoff C, Karim R, Ahmed E, Camara O, Rhode KS.
> J Interv Card Electrophysiol. 2017 Sep 7. doi: 10.1007/s10840-017-0281-3.

> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5633640/

Please cite this reference when using this code.

## Pipeline

The pipeline is split in four parts.

* run_standardization: standardizes meshes from a raw mesh. Depends on [VMTK], [VTK] and [MeshLab]
* run_currents: registers mesh to an atlas using currents registration. Depends on [VTK], [MATLAB] and ./currents_build
* run_sum: computes standardized unfold map. Depends on [VMTK] and [VTK]
* run_quantification: computes regional quantification (extent or mean value). Depends on [VMTK] and [VTK]


## Instructions

Clone the repository and cd into it:
```sh
$ git clone https://github.com/catactg/sum
$ cd sum
```

Set paths according to your system in `constants.py`
```
MATLAB_BIN_PATH
CURRENTS_BUILD_PATH
MESHLABSERVER_PATH
```
Default parameters of the algorithm can be modified in `main.py` inside the `get_parameters()` method.

## Usage

```
main.py [-h] --meshfile MESHFILE --datatype DATATYPE
               [--pvcliptype PVCLIPTYPE]
               [--paired_unfold_disk PAIRED_UNFOLD_DISK]
               [--mitral_clip_type MITRAL_CLIP_TYPE] [--use_seed_selector]
               [--use_laa_seed] [--use_similarity] [--use_glyphs]
               [--visualize] [--skip_standardization] [--skip_currents]
               [--skip_sum] [--skip_quantification]


  -h, --help              show this help message and exit
  --meshfile MESHFILE     Full path to mesh file in VTK format
  --datatype DATATYPE     Data type to process: force | lge | lat
  --pvcliptype PVCLIPTYPE How to clip the pulmonary veins: short | long
  --paired_unfold_disk PAIRED_UNFOLD_DISK
                          Full path to a paired unfold disk. The edges of the
                          current unfold will be overlaid on the paired unfold
                          disk
  --mitral_clip_type MITRAL_CLIP_TYPE
                          The algorithm will compute a mitral plane based on the
                          body and PVs centroids. If the meshfile already
                          contains a mitral plane clip, use the 'manual' option
                          to preserve it
  --use_seed_selector     Seed selection will always run for a new case. To
                          select new seeds, activate this flag
  --use_laa_seed          Activate seed for appendage
  --use_similarity        The method uses an affine transform to initialize mesh
                          registration. Activate this flag to use a similarity
                          transform
  --use_glyphs            Activate glyph on visualization. Recommended for spare
                          measurements (e.g. some lat meshes)
  --visualize             Activate visualization
  --skip_standardization  Skip run_standardization step
  --skip_currents         Skip run_currents step
  --skip_sum              Skip run_sum step
  --skip_quantification   Skip run_quantification step
```

## Test data
Data was created by using a mesh from the [LASC] challenge and adding mock scalars onto it: force, lge, and lat.

Run the command lines below and the PNG output should be identical to the one in `./data/expected_output`.

```
python ./main.py --meshfile ./data/force/mock_force.vtk --datatype force

python ./main.py --meshfile ./data/lge/mock_lge.vtk --datatype lge --paired_unfold_disk ./data/force/mock_force_FTI_disk_uniform.vtp

python ./main.py --meshfile ./data/lat/mock_lat.vtk --datatype lat

python ./main.py --meshfile ./data/lat/mock_lat.vtk --datatype lat --use_glyph --skip_standardization --skip_currents
```

## Dependencies

The scripts in this repository were successfully run with:
- [Python] 2.7.14
- [NumPy] 1.13.3
- [VMTK] 1.4
- [VTK] 8.0.1
- [MATLAB] R2017b
- [MeshLab] 1.3

### Python packages installation

The easiest way to install everything you need using the VMTK [conda] package. It inclues Python, NumPy, VMTK and VTK.

First, download miniconda and install it from terminal:
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
sh Miniconda2-latest-MacOSX-x86_64.sh
```

Create an environment and activate it:
```
conda create --name vmtk
source activate vmtk

```

Finally, install vmtk:
```
conda install -c vmtk/label/dev itk vtk vmtk
```

Anytime you wish to use vmtk (e.g. to run the code in this repository) you will need to activate the environment:
```
source activate vmtk
```

[LASC]:http://github.com/catactg/lasc
[Python]:http://www.python.org
[NumPy]:http://www.numpy.org
[VMTK]:http://www.vmtk.org
[VTK]:http://www.vtk.org
[MATLAB]:http://www.mathworks.com
[MeshLab]:http://www.meshlab.net
[conda]:https://conda.io/docs/


## License

BSD 2-Clause