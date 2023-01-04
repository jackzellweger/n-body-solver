Function Declarations

```
PythagoreanTriple[distance_] := Module[{x, y, a, q, r},
  x = RandomReal[{0, distance}];
  a = Quiet[Solve[distance == Sqrt[x^2 + y^2], y, Reals]];
  {q, r} = {x, Values @@ a[[2]]};
  {q, r} = RandomChoice[{{q, r}, {q, -r}, {-q, r}, {-q, -r}}];
  Return[{q, r}]
  ]
```

# The N-Body Solver

### Introduction

I've always been fascinated by the n-body problem. The way that objects dance around one-another is beautiful. Their chaos, and the way their movements hang in such delicate balance is fascinating.

But amazingly, all these beautifully complex interactions just follow one simple rule:

$$ F = G\frac {m} {a} $$

That's right, like 99% of physics, it all comes from Newton's Second Law. Almost everything in our universe moves—on a macro scale—according to this law. These movements depend on such a simple equation; it should be easy to write a program that calculates them, right?

Since I have a physics background, let's draw inspiration from the scale and movements of our own solar system. I will use numbers that roughly approximate (to +/- 3 orders of magnitude) the sizes and masses of the planets.

There are a few constants, then, that we need to have handy for our program to access. We will use combinations of these constants to define other variables that we will use to run our program. For example, we will use a combination of constants to define the orbital period of the Earth, T. We will then use T to define the step size for our program. Everything has to agree on the same order of magnitude for this program to work—everything has to be in-sync.

**Note: All units are given in SI units.**

### Declaring Celestial Variables & Constants

The Universal Gravitational Constant (UGC). This constant defines the strength of gravity per unit mass, relating kilograms and newtons:

`G = 6.67 10^-11;`

The mass of the sun in kilograms:

`ms = 1.988435 10^30;`

The mass of the earth in kilograms:

`me = 5.97219*10^24;`

The radius of Sun in meters:

`R = 6.955*10^8;`

Earth's starting x- coordinate in meters:

`xn = 1.506*^11;`

Earth's starting y-coordinate in meters:

`yn = 0.;`

Earth's starting x-velocity in meters per second:

`vxn = 0.;`

Starting speed of earth in y direction in meters per second

`vyn = 29598.;`

The magnitude of the semi-major axis in meters

`a = (G ms*(xn))/(2 G ms - xn vyn^2);`

Here, we use Kepler's third law, and some of the constants we've already defined to calculate the orbital period of Earth. Why do we need that? It's just nice to have; it gives us a nice to-scale unit we can use to step-through our results.

```T = (2 a^(3/2) (me + ms) \[Pi])/(Sqrt[G] ms^(3/2));```

This is the approximate distance from Earth to the Sun at the time the program is run:

```
earthDistance = 
  QuantityMagnitude[
   UnitConvert[
    Entity["Star", "Sun"][
     EntityProperty["Star", "DistanceFromEarth", {"Date" -> Now}]], 
    "Meters"]];
```

The Earth's Angular Velocity in radians per second:

```
earthAngularVel = 
  QuantityMagnitude[
   UnitConvert[
    2 \[Pi]/Entity["Planet", "Earth"][
       EntityProperty["Planet", "OrbitPeriod"]]]];
```

We can calculate the Earth's speed around the Sun from its orbital distance and orbital angular velocity:

`earthSIVel = earthDistance*earthAngularVel;`

The Earth's mass in kilograms:

```
me = QuantityMagnitude[
   UnitConvert[Entity["Planet", "Earth"][EntityProperty["Planet", "Mass"]]]];
```

The Sun's mass in kilograms:

```
ms = QuantityMagnitude[
   UnitConvert[Entity["Star", "Sun"][EntityProperty["Star", "Mass"]]]];
```

### User Input

Change this number to the number of planets we'd like to simulate in the system:

`planetNumber = 3;`

The distances of the bodies from the Sun in astronomical units (AU).
- Change these values to whatever you want.
- The number of entries here has to be equal to `planetNumber`.

```
planetDistances = {5, 5, 5, 5, 5};
If[planetNumber != Length[planetDistances], Abort[]];
```

Here, we create an array to hold each planet's distance in meters:

```
distanceArrayAU = planetDistances*earthDistance;
If[planetNumber != Length[distanceArrayAU], Abort[]];
```

We set the mass of each of the planets to a random value between 100,000x and 111,000x the Earth's mass; we then put them in an array:

```
massArray1 = RandomReal[{100000, 110000}, planetNumber];
Assert[Length[massArray1] == Length[planetNumber]];
If[Length[massArray1] != planetNumber, Abort[]]
earthMassesArray = massArray1*me;
```

A more descriptive variable name would be nice, but I chose `m` as the name here because it shortens the equation that we eventually have to plug this into.

`m = earthMassesArray;`

This is where we used the function we declared in the function declaration section. Here, we solve the Pythagorean Theorem to find our starting X and Y values on the Cartesian plane. In the equation a^2+b^2=c^2, we set the distance the user has chosen equal to c (the hypotenuse), and solve for a & b using a random seed.

```
coordsTable = Table[PythagoreanTriple[earthDistance], 3];
initialX = coordsTable[[All, 1]];
initialY = coordsTable[[All, 2]];
```

Then we do the same thing with velocity vectors, this time using Earth's velocity as c, the hypotenuse, and finding two variables  a and b that work.

`coordsTable = Table[PythagoreanTriple[earthSIVel], 3];`

Here, we prime our position and velocity variables for use in the equations of motion.

```
xarr = Table[x[n][t], {n, 1, planetNumber}];
yarr = Table[y[n][t], {n, 1, planetNumber}];
vxarr = Table[vx[n][t], {n, 1, planetNumber}];
vyarr = Table[vy[n][t], {n, 1, planetNumber}];
```

And here we calculate velocities:

```
vxarr = coordsTable[[All, 1]]/3;
vyarr = coordsTable[[All, 2]]/3;
```

Now for a bit of theory. The gravitational force between n planets is being calculated using Newton's law of universal gravitation.

$$ \vec {F} = G * \frac {m_ {1} * m_ {2}} {r^{2}} * \hat {r} $$

 We can project this vector into its two one-dimensional components\[LongDash]one for the x-dimension and one for the y-dimension, in order to calculate the force in each direction independently. We get....

$$ x_{j}''(t) = \sum_{i=1}^{\text{planetNumber}} (1-\delta_{ij}) \frac{G * m_{j} * (x_{i}(t) - x_{j}(t))}{\sqrt{(x_{i}(t) - x_{j}(t))^2 + (y_{i}(t) - y_{j}(t))^2}} $$


and

$$ (y_ {j}'')(t) = \sum_{i=1}^{\text{planetNumber}} (1-\delta_{ij}) \frac{G * m_{j} * (y_{i}(t) - y_{j}(t))}{\sqrt{(x_{i}(t) - x_{j}(t))^2 + (y_{i}(t) - y_{j}(t))^2}} $$

If i and j are the same, the force is set to 0 (since a planet doesn't exert a gravitational force on itself). If i and j are different, the force is calculated using the equation for gravitational force. We use the Kronecker Delta to represent this above.

Here's the above equations converted into a form that Mathematica can read.

```
equationsOfMotionX = Table[(x[j]'')[t] == Sum[If[i == j, 0, G*m[[j]]*((x[i])[t] - (x[j])[t])/Sqrt[((x[i])[t] - (x[j])[t])^2 + ((y[i])[t] - (y[j])[t])^2]*(((x[i])[t] - (x[j])[t])^2 + ((y[i])[t] - (y[j])[t])^2)], {i, 1, planetNumber}], {j, 1, planetNumber}];

```
Since this is really dense, this is what they look like in LaTex

$$
\begin{aligned}
equationsOfMotionX = \left\{ \left( x_j'' \right)(t) = \sum_{i=1}^{planetNumber} \left[ \begin{aligned} 0 &\quad\text{if }i=j\\ G\ m_j \frac{\left( x_i \right)(t) - \left( x_j \right)(t)}{\sqrt{\left( x_i \right)(t) - \left( x_j \right)(t)^2 + \left( y_i \right)(t) - \left( y_j \right)(t)^2}} \left( \left( x_i \right)(t) - \left( x_j \right)(t)^2 + \left( y_i \right)(t) - \left( y_j \right)(t)^2 \right) \right] \right\} \quad\text{for }j\in\left\{ 1,\dots,planetNumber \right\}
\end{aligned}
$$














