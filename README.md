# The N-Body Solver

## Introduction

I've always been fascinated by the n-body problem. The way that objects dance around one-another is beautiful. Their chaos, and the way their movements hang in such delicate balance is fascinating.

But amazingly, all these beautifully complex interactions just follow one simple rule:

$$ F = m a $$

That's right, like 99% of physics, it all comes from Newton's Second Law. Almost everything in our universe moves—on a macro scale—according to this law. These movements depend on such a simple equation; it should be easy to write a program that calculates them, right?

Since I have a physics background, let's draw inspiration from the scale and movements of our own solar system. I will use numbers that roughly approximate (to +/- 3 orders of magnitude) the sizes and masses of the planets.

There are a few constants, then, that we need to have handy for our program to access. We will use combinations of these constants to define other variables that we will use to run our program. For example, we will use a combination of constants to define the orbital period of the Earth, `T`. We will then use `T` to define the step size for our program. Everything has to agree on the same order of magnitude for this program to work—everything has to be in-sync.

**Note: All units are given in SI units.**


## Declaring Celestial Variables & Constants

The Universal Gravitational Constant (UGC). This constant defines the strength of gravity per unit mass, relating kilograms and newtons:

`G = 6.67 10^-11;`

The mass of the Sun in kilograms:

`ms = 1.988435 10^30;`

The mass of the Earth in kilograms:

`me = 5.97219*10^24;`

The radius of Sun in meters:

`R = 6.955*10^8;`

Earth's starting `x`-coordinate in meters:

`xn = 1.506*^11;`

Earth's starting `y`-coordinate in meters:

`yn = 0.;`

Earth's starting `x`-velocity in meters per second:

`vxn = 0.;`

Starting speed of earth in `y`-direction in meters per second

`vyn = 29598.;`

The magnitude of the semi-major axis in meters

`a = (G ms*(xn))/(2 G ms - xn vyn^2);`

Here, we use Kepler's third law, and some of the constants we've already defined to calculate the orbital period of Earth. Why do we need that? It's just nice to have; it gives us a nice to-scale unit we can use to step-through our results.

```T = (2 a^(3/2) (me + ms) \[Pi])/(Sqrt[G] ms^(3/2));```

This is the approximate distance from Earth to the Sun at the time the program is run:

```mathematica
earthDistance = 
  QuantityMagnitude[
   UnitConvert[
    Entity["Star", "Sun"][
     EntityProperty["Star", "DistanceFromEarth", {"Date" -> Now}]], 
    "Meters"]];
```

The Earth's Angular Velocity in radians per second:

```mathematica
earthAngularVel = 
  QuantityMagnitude[
   UnitConvert[
    2 \[Pi]/Entity["Planet", "Earth"][
       EntityProperty["Planet", "OrbitPeriod"]]]];
```

We can calculate the Earth's speed around the Sun from its orbital distance and orbital angular velocity:

`earthSIVel = earthDistance*earthAngularVel;`

The Earth's mass in kilograms:

```mathematica
me = QuantityMagnitude[
   UnitConvert[Entity["Planet", "Earth"][EntityProperty["Planet", "Mass"]]]];
```

The Sun's mass in kilograms:

```mathematica
ms = QuantityMagnitude[
   UnitConvert[Entity["Star", "Sun"][EntityProperty["Star", "Mass"]]]];
```

## Declaring Functions

This is a function I wrote that will prove helpful to set intitial conditions.

```mathematica
PythagoreanTriple[distance_] := Module[{x, y, a, q, r},
  x = RandomReal[{0, distance}];
  a = Quiet[Solve[distance == Sqrt[x^2 + y^2], y, Reals]];
  {q, r} = {x, Values @@ a[[2]]};
  {q, r} = RandomChoice[{{q, r}, {q, -r}, {-q, r}, {-q, -r}}];
  Return[{q, r}]
  ]
```

## Getting User Input

Change this number to the number of planets we'd like to simulate in the system:

`planetNumber = 3;`

The distances of the bodies from the Sun in astronomical units (AU).
- Change these values to whatever you want.
- The number of entries here has to be equal to `planetNumber`.

```mathematica
planetDistances = {5, 5, 5, 5, 5};
If[planetNumber != Length[planetDistances], Abort[]];
```

## The Data Structures

Here, we create an array to hold each planet's distance in meters:

```mathematica
distanceArrayAU = planetDistances*earthDistance;
If[planetNumber != Length[distanceArrayAU], Abort[]];
```

We set the mass of each of the planets to a random value between 100,000x and 111,000x the Earth's mass; we then put them in an array:

```mathematica
massArray1 = RandomReal[{100000, 110000}, planetNumber];
Assert[Length[massArray1] == Length[planetNumber]];
If[Length[massArray1] != planetNumber, Abort[]]
earthMassesArray = massArray1*me;
```

A more descriptive variable name would be nice, but I chose `m` as the name here because it shortens the equation that we eventually have to plug this into.

`m = earthMassesArray;`

This is where we used the function we declared in the function declaration section. Here, we solve the Pythagorean Theorem to find our starting `X` and `Y` values on the Cartesian plane. In the equation $a^2+b^2=c^2$, we set the distance the user has chosen equal to `c` (the hypotenuse), and solve for `a` and `b` using a random seed.

```mathematica
coordsTable = Table[PythagoreanTriple[earthDistance], 3];
initialX = coordsTable[[All, 1]];
initialY = coordsTable[[All, 2]];
```

Then we do the same thing with velocity vectors, this time using Earth's velocity as `c`, the hypotenuse, and finding two variables `a` and `b` that work.

`coordsTable = Table[PythagoreanTriple[earthSIVel], 3];`

Here, we prime our position and velocity variables for use in the equations of motion.

```mathematica
xarr = Table[x[n][t], {n, 1, planetNumber}];
yarr = Table[y[n][t], {n, 1, planetNumber}];
vxarr = Table[vx[n][t], {n, 1, planetNumber}];
vyarr = Table[vy[n][t], {n, 1, planetNumber}];
```

And here we calculate velocities:

```mathematica
vxarr = coordsTable[[All, 1]]/3;
vyarr = coordsTable[[All, 2]]/3;
```

## Calculating Forces Between Bodies

Now for a bit of theory. The gravitational force between `n` planets is being calculated using Newton's Law of Gravity.

$$ \vec{F} = G * \frac{m_{1} * m_{2}} {r^{2}} * \hat{r} $$

 We can generalize this vector equation for `n` bodies, and project it its two one-dimensional components—one for the `x`-dimension and one for the `y`-dimension in order to calculate the force in each direction independently.

$$ x_{j}''(t) = \sum_{i=1}^{\text{planetNumber}} (1-\delta_{ij}) \frac{G * m_{j} * (x_{i}(t) - x_{j}(t))}{\sqrt{(x_{i}(t) - x_{j}(t))^2 + (y_{i}(t) - y_{j}(t))^2}} $$

and

$$ y_ {j}''(t) = \sum_{i=1}^{\text{planetNumber}} (1-\delta_{ij}) \frac{G * m_{j} * (y_{i}(t) - y_{j}(t))}{\sqrt{(x_{i}(t) - x_{j}(t))^2 + (y_{i}(t) - y_{j}(t))^2}} $$

If `i` and `j` are the same, the force is set to 0 (since a planet doesn't exert a gravitational force on itself). If `i` and `j` are different, the force is calculated using the equation for gravitational force. We use the Kronecker Delta to represent this above.

Let's convert the above equations into a form that *Mathematica* can read. Here is the equation of motion for the `x` dimension:

```mathematica
equationsOfMotionX = 
  Table[
    (x[j]'')[t] == 
      Sum[
        If[i == j, 0, 
          G*m[[j]]*
          ((x[i])[t] - (x[j])[t])/
          Sqrt[
            ((x[i])[t] - (x[j])[t])^2 + 
            ((y[i])[t] - (y[j])[t])^2]*
          (((x[i])[t] - (x[j])[t])^2 + 
           ((y[i])[t] - (y[j])[t])^2)
        ], 
        {i, 1, planetNumber}
      ], 
      {j, 1, planetNumber}
    ];
```

Here is the equation of motion for the `y` dimension:

```mathematica
equationsOfMotionY = 
  Table[
    (y[j]'')[t] == 
      Sum[
        If[i == j, 0, 
          G*m[[j]]*
          ((y[i])[t] - (y[j])[t])/
          Sqrt[
            ((x[i])[t] - (x[j])[t])^2 + 
            ((y[i])[t] - (y[j])[t])^2]*
          (((x[i])[t] - (x[j])[t])^2 + 
           ((y[i])[t] - (y[j])[t])^2)
        ], 
        {i, 1, planetNumber}
      ], 
      {j, 1, planetNumber}
    ];
```

This code defines two sets of equations: one for the `x`-coordinates of the planets (`equationsOfMotionX`) and one for the `y`-coordinates (`equationsOfMotionY`).

- `x[j][t]` and `y[j][t]` are the variables that representing the `x` & `y` coordinates of the `j`th planet at time `t`.
- `x[i][t]` and `y[i][t]` are variables representing the `x` and `y` coordinates, respectively, of the `i`th planet at time `t`.
- `G` is the gravitational constant.
- `m[j]` is the mass of the `j`th planet.

In place of the Kronecker Delta, we use Mathematica's `If` function to check whether the planet being considered (`i`) is the same as the planet whose equation of motion is being defined (`j`). If it is, we don't include that force since planets don't exert forces on themselves.

## Initial Conditions

We combine the equations of motion into a big list:

```mathematica
equationsOfMotion = Join[equationsOfMotionX, equationsOfMotionY]
```

We then use the `Table` function to define a list of initial conditions for the positions and velocities of the planets at time `t = 0`. Notice the 0s plugged into the functions. The initial positions are in `initialXArr` and `initialYArr`, while the initial velocities are in `initialXDerivArr` and `initialYDerivArr`.

```mathematica
initialXArr = Table[x[j][0], {j, 1, planetNumber}];
initialYArr = Table[y[j][0], {j, 1, planetNumber}];
initialXDerivArr = Table[Derivative[1][x[j]][0], {j, 1, planetNumber}];
initialYDerivArr = 
  Table[Derivative[1][y[j]][0], {j, 1, planetNumber}];
```

We use`Thread` function to create a list of equations for the numerical solver by pairing  of the initial position and velocity lists with the corresponding functions. For example, `Thread[xarr==initialX]` results in `{y[1][t]==6.96*10^10,y[2][t]==-1.5*10^11,y[3][t]==5.6*10^10}`. This equation is now ready for us to substitute 0 for `t` in all cases.


So, we substitute `t = 0` into the resulting list of equations. We then store the resulting list of equations, which represents the system of differential equations to be solved. We then do a little list clean-up.

```mathematica
initialConditions = 
  Flatten[Join[{Thread[xarr == initialX], Thread[yarr == initialY], 
     Thread[initialXDerivArr == vxarr], 
     Thread[initialYDerivArr == vxarr]}]];
equationsForNDSolve = 
  Join[equationsOfMotion, initialConditions /. {t -> 0}];
allVars = Flatten[{xarr, yarr}];
```

## Numerically Solving The Equation of Motion

The code then uses Mathematica's `NDSolve` function to numerically solve the system of differential equations defined by `equationsOfMotion` and the initial conditions, over a time interval from `t = 0` to `t = 25 * T`, where `T` is the period of the outermost planet.


```mathematica
allOrbits = 
   NDSolve[equationsForNDSolve, allVars, {t, 0, 25*T}(*, 
    AccuracyGoal \[Rule] 20, 
    Method\[Rule]"ExplicitRungeKutta"*)]; // Quiet
```

This generates the list of equations in the right order to plug in to our `ParametricPlot` function.

```mathematica
drawList = {Nothing};
For[i = 1, i <= planetNumber, i++,
  AppendTo[drawList, xarr[[i]]];
  AppendTo[drawList, yarr[[i]]];
  ];
```

Then we the paths of the planets parametrically. This long line of code separates the plotting (`ParametricPlot`) from the plot image-rasterization (`Rasterize`) step. In other words, for each frame of the orbit calculation, we plot the path parametrically, then rasterize it. Then we use `ListAnimate` to step through the rasterize images. This is a common technique to increase the performance of animations.

```mathematica
images = Table[
   Rasterize[
    ParametricPlot[drawList /. allOrbits, {t, 0, n*T}, 
     AspectRatio -> 1., ImageSize -> Large, 
     PlotRange -> {{-5*10^11, 5*10^11}, {-5*10^11, 5*10^11}}, 
     ColorFunction -> (ColorData["Rainbow"][#3] &), Axes -> False], 
    ImageResolution -> 150], {n, 0.01, 5, 0.03}];
ListAnimate[images, AnimationRunning -> True, AnimationRate -> 24, 
 ImageSize -> 450]
```

We then export the animation of the orbiting bodies...

`Export["/path/to/file.mov", images]`


## Result

https://user-images.githubusercontent.com/20098240/210496445-fa02860f-a783-4984-8f88-cf02956dd801.mov
