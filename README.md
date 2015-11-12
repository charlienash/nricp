# nricp - Non-rigid iterative closest point

https://github.com/charlienash/nricp

nricp is a MATLAB implementation of a non-rigid variant of the iterative closest point algorithm. It can be used to register 3D surfaces or point-clouds. The algorithm is an implementation of the method described in the following paper along with some additional features:

'Optimal Step Nonrigid ICP Algorithms for Surface Registration', Amberg, Romandhani and Vetter, CVPR, 2007.

Features:
* Non-rigid and local deformations of a template surface or point cloud.
* Iterative stiffness reduction allows for global intitial transformations that become increasingly localised.  
* Optional initial rigid registration using standard iterative closest point.
* Optional bi-directional distance metric which encourages surface deformations to cover more of the target surface

Limitations:
* Does not yet support missing data in either the target or source surfaces.

## Dependencies

Requires:
* geom3d - http://uk.mathworks.com/matlabcentral/fileexchange/24484-geom3d
* Toolbox Graph - http://uk.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph
* Iterative Closest Point - http://uk.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point 

## Installation

Download the nricp directory and the dependencies and add them to your MATLAB path. 

## Attribution

If you use this implementation in your academic projects, please cite the paper by Amberg et al:

```bibtex
@inproceedings{amberg2007optimal,
  title={Optimal step nonrigid icp algorithms for surface registration},
  author={Amberg, Brian and Romdhani, Sami and Vetter, Thomas},
  booktitle={Computer Vision and Pattern Recognition, 2007. CVPR'07. IEEE Conference on},
  pages={1--8},
  year={2007},
  organization={IEEE}
}
```

## Contact
charlie.tc.nash@gmail.com

