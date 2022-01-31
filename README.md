# Materials Science Animations
A home for interesting and intuitive animations related to materials science.

## Example Gallery
### Miscibility gap from the common tangent construction
This animation shows how we can get a miscibility gap by constructing a common
tangent line and tracking the corresponding compositions as the temperature increases.
![](./media/videos/miscibility_gap/720p30/MiscibilityGap_ManimCE_v0.14.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/thermodynamics/miscibility_gap.py#L123)


### Fourier Series as a solution to the heat equation
This animation explains how we can approximate a square wave by incrementally
adding up sine curves.
![](./media/videos/heat_equation/720p30/FourierSeries_ManimCE_v0.14.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/kinetics/heat_equation.py#L5)

### Three-fold rotational symmetry of a cube
This animation shows the three-fold rotational symmetry of a cube.
![](./media/videos/cube_symmetry/720p30/CubeSymmetry_ManimCE_v0.14.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/crystallography/cube_symmetry.py#L4)


### Calculate the time cone volume in 1D, 2D, 3D space-time
An animation that visually explains the calculation steps 
for the time cone volume in 1D, 2D, 3D space-time.
![](./media/videos/present/720p30/CalcTimeConeVolume_ManimCE_v0.12.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/present.py#L1293)

### Phase fraction to Time-Temperature-Transformation diagram
An animation that shows how we can get a Time-Temperature-Transformation phase diagram 
by plotting the (volume) phase fraction transformed as a function of time at different temperatures.
![](./media/videos/present/720p30/PhaseFractionToTTT_ManimCE_v0.12.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/present.py#L1709)

### Convert a probability to an overall phase fraction
![](./media/videos/present/720p30/ProbToPhaseFraction_ManimCE_v0.12.0.gif)
[[Source code here]](https://github.com/rpw199912j/matsci_animation/blob/master/present.py#L859)

This animation shows, at an arbitrary time, how we can get the overall phase fraction transformed/untransformed by sampling the
probability of being transformed/untransformed at multiple points in space. 
Here I assume the material undergoing the phase transformation is nucleating and growing on a 2D plane. 
Since I set the probability of being untransformed at an arbitrary point in space 
(or the probability that no nucleation has occurred) to be 0.2, the overall phase fraction transformed should be 0.8.
What I first showed is that we are sampling evenly a 10 by 10 grid on a 2D plane, giving us 100 squares. 
Within each square on that grid, there is a 0.2 probability that it will not be transformed 
(i.e., 0.8 probability that it will be transformed). 
If a phase transformation does occur within that square, 
we populate that square with a white circle; if there isn't a phase transformation, 
that square will be left blank. 
To simulate the phase transformation process, 
I use a random number generator to give me a 10 by 10 grid of numbers in the range of [0, 1). 
For each number, if it is greater than or equal to 0.2, then a phase transformation has occurred; 
if the number is less than 0.2, then a phase transformation has not occurred. 
After all of this, the number of white circles should be close to 80. 
However, since we only have a 10 by 10 grid, and also due to the probabilistic nature of the simulation, 
the number of white circles won't always be exactly 80. 
We could, in theory, increase the size of the sampling grid and let it approach a number on the order of 10<sup>23</sup>, 
which would mimic a real-life material and also give us a fraction very close to 0.8. 
The last part of the animation is just fun to visualize. 
What I am trying to show is that instead of counting the number of white circles and divide it by 100 to get 
the overall phase fraction transformed, we can treat all the circles as rigid bodies and simulate the effect of gravity 
so that all the circles should fill up to about the 80% horizontal line. 