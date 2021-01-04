# cop-sim-demo
Some examples of contact simulation using the COP resolution method.

##### 01-single-contact:
Simulation of a falling and rolling capsule on a plane surface.
Note: This does not use the cop_simulator library. Instead, it directly resolves the
COP point for a line contact using the `resolveCOPLineContactWithLastCOPSol` function.
For collision resolution, the `solveCollLCPOnePoint` and `solveCollLCPPoint` functions
are used.

##### 02-two-capsules:
Simulation of several falling, colliding and rolling capsules using the cop_sim framework.
Note: This simulation fails currently if a simultaneous contact occurs between one capsule
with the ground as well as with another cylinder lying flat on the ground. i.e. only upto
1 simultaneous point contact and 1 line contact are resolved.
Note: 02-two-capsules/collisionDetectionDemo.cpp: This visualizes the minimum distance
computation between two capsules with a single collision between them.

##### 03-one-cylinder:
Simulation of a surface contact between one cylinder and a ground surface.
A time-increasing force is applied to the top surface of the cylinder causing it to topple
over due to friction at the surface contact.

##### 04-two-arm-rolling:
Simulation of a robot arm rolling a rolling pin with the palm of its hand. The contact surface of the palm is modeled as a single capsule perpendicular to the rolling pin. As a result, there is only a simultaneous 1 point- 1 line contact. Eventually, this will be extended to two arms with 3 simultaneous line contacts.

##### 05-one-box:
Simulation of a surface contact between one box and a ground surface.

##### 06-one-box:
Simulation of a surface contact between one n-sided regular pyramid and a ground surface.
