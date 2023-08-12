# cco_vasculature

## WIP Implementation of vasculature tree generation based on the 1993 paper of Schreiner et al.

### Model Constants(constants.py)
The model constants can be changed as needed.  
###### Description

SEGMENT\_ASPECT\_RATIO - The minimum ratio of the length of the vessel segment to its radius

WALL\_THICKNESS\_LUMEN\_RATIO - Ratio of wall thickness to the segment radius

OPENSCAD\_SCALING\_FACTOR - Only applies to exporter, scales the model so that it is not too small

CORE\_COUNT - Constant for WIP parallelism mechanism. Parallel processing is currently not implemented.


### Usage
- Set constants in constants.py
- Run 'cli.py'
- Exported output is stored in 'test.scad'
- Open 'test.scad' in OpenSCAD


