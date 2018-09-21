# 3D DENSE Plugin for Crescent Organ
[*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) affiliate plugin for analyzing 3D Cine **d**isplacement **en**coding with **s**timulated **e**choes (DENSE) MRI data of a crescent-shaped organ

## Applications
- Right ventricle of small animals like rodent, which have a much more crescent-shaped cavity
- The heart in which the biventricular shape is sharper at the right ventricle insertion points

## Features
* Iteration Algorithm for Adjusting the Finite Element Mesh into its Actual Crescent Shape
* Insertion Tissue Resection
* A Graphical User Interface for User-defined Noise Suppression
* A Graphical User Interface for Visualizing the Deforamtion of Finite Elements of Interest

## Requirements
- [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) 
	- The core software
- The biventricular branch of [*DENSE3D Plugin*](https://github.com/suever/dense3D_plugin/tree/biventricular)
	- Call for the scripts under the private directory inside it

## Compatibility
Compatible with [v0.5.0](https://github.com/denseanalysis/denseanalysis/blob/master/CHANGELOG) of [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis)

## Installation
After installation of [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) and the biventricular branch of [*DENSE3D Plugin*](https://github.com/suever/dense3D_plugin/tree/biventricular), run the following from the MATLAB command line:

```matlab
plugins.PluginManager.import('https://github.com/MMoTH/DENSE3D_Plugin_4CrescentOrgan')
```

## Credits
* This package was created with the [*denseanalysis_plugin_demo*](https://github.com/denseanalysis/denseanalysis_plugin_demo) project template.
* Matlab scripts [`wrapMesh.m`], [`circumferentialParameterize.m`],  and [`computeStrains.m`] are modified from the [*DENSE3D Plugin*](https://github.com/suever/dense3D_plugin/tree/biventricular) project.

## Known Bugs

## Known Issues
1.Error in plugins.PluginMenu/checkAvailability

Temperal Solution: Disable checkAvailability
