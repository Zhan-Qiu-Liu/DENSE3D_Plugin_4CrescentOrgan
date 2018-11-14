## 3D DENSE Plugin for Crescent Organ
[*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) affiliate plugin for analyzing 3D Cine **d**isplacement **en**coding with **s**timulated **e**choes (DENSE) MRI data of a crescent-shaped organ

### Applications
- Right ventricle of small animals like rodent, which have a much more crescent-shaped cavity
- The heart in which the biventricular shape is sharper at the right ventricle insertion points

### Features
* Iteration Method for Adjusting the Finite Element Mesh into its Actual Crescent Shape
* Precise Segmentation of Right Ventricular Free Wall from Left Ventricular Free Wall (by taking the fact into account that the anterior and inferior insertions of the RV to the LV shift circumferentially along the ventricular long axis)
* A Graphical User Interface for User-defined Noise Suppression
* A Graphical User Interface for Selecting the Finite Elements of Interest and Visualizing their Deforamtion 

### Requirements
- [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) 
	- The core software for 2D contouring
- The biventricular branch of [*DENSE3D Plugin*](https://github.com/suever/dense3D_plugin/tree/biventricular)
	- Call for the scripts under the private directory inside it

### Compatibility
Compatible with [v0.5.0](https://github.com/denseanalysis/denseanalysis/blob/master/CHANGELOG) of [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis)

### Installation
After installation of [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) and the biventricular branch of [*DENSE3D Plugin*](https://github.com/suever/dense3D_plugin/tree/biventricular), run the following from the MATLAB command line:

```matlab
plugins.PluginManager.import('https://github.com/MMoTH/DENSE3D_Plugin_4CrescentOrgan')
```

>**NOTE**: For analysis only, click once `Initialize DENSE3DPlugin4CrescentOrgan` under the menu tab `Plugin_DENSE3D4CrescentOrgan` every time [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) is ran.

### Credits
- This package was created with the [*denseanalysis_plugin_demo*](https://github.com/denseanalysis/denseanalysis_plugin_demo) project template.
- Copyright (c) of the following improved Matlab scripts belongs to 2016 [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) Contributors and Zhanqiu Liu <lafeir.lew@gmail.com>:
	- [`wrapMesh.m`](wrapMesh.m), [`circumferentialParameterize.m`](circumferentialParameterize.m),  [`computeStrains.m`](computeStrains.m), [`analysisFcnREPL.m`](analysisFcnREPL.m), [`spl2strainFcnREPL.m`](spl2strainFcnREPL.m), [`exportMatREPL.m`](exportMatREPL.m), [`MeshControl.m`](MeshControl.m), and [`linearInterp.m`](linearInterp.m)
- Copyright (c) of the following orignial Matlab scripts soley belongs to 2016 [*DENSEanalysis*](https://github.com/denseanalysis/denseanalysis) Contributors:
	- [`contourDataREPL.m`](contourDataREPL.m), [`computeTransmuralDataREPL.m`](computeTransmuralDataREPL.m), [`resetFcnREPL.m`](resetFcnREPL.m), [`switchstateREPL.m`](switchstateREPL.m), and [`animatedgif.m`](animatedgif.m)

### Attribution

You are welcome to use the *DENSE3D_Plugin_4CrescentOrgan* in your research provided that you include the following citation in your work:

>Zhanqiu (2018). DENSE3D_Plugin_4CrescentOrgan. Retrieved
>from https://github.com/MMoTH/DENSE3D_Plugin_4CrescentOrgan


### Known Bugs

### Known Issues
1.Error in plugins.PluginMenu/checkAvailability

Temperal Solution: Disable checkAvailability
