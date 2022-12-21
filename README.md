# Galaxy-Collisions
MATLAB Galaxy Collision Simulator

The aim of this project is to use MATLAB code and finite difference approximations (FDAs) to implement a simplified model of galactic collisions. The exercise was inspired by the Alan Toomre model, in which each galaxy is modelled as a central particle (known as a core) with some gravitating mass m. A core is orbited circularly by n stars, which conversely have vanishing mass — i.e. they experience acceleration due to gravity from the core(s), but they themselves do not induce gravitational acceleration in any other body.

Using the finite difference model, we take a function on a continuous time domain 0 ≤ t ≤ tmax, and approximate it by converting to a discrete time domain. The number of time steps, as well as the size of the time steps in our discrete domain, are determined by the chosen level parameter, ℓ:

<img width="201" alt="FDA1" src="https://user-images.githubusercontent.com/48260444/209003598-576be52f-1da0-404c-8e54-28c938c9a12f.png">
<img width="372" alt="FDA2" src="https://user-images.githubusercontent.com/48260444/209003616-f6144972-dc3e-498b-a33b-6d6499fe0a8a.png">
<img width="529" alt="FDA3" src="https://user-images.githubusercontent.com/48260444/209003629-bbda030c-49fe-447d-a101-87f785bc303e.png">

Below are the essential equations for this model:

<img width="303" alt="EQN1" src="https://user-images.githubusercontent.com/48260444/209003256-11254202-bcbc-4ae1-b269-32a90a0f3959.png">
<img width="468" alt="EQN2" src="https://user-images.githubusercontent.com/48260444/209003363-466206c1-f0d0-43c7-b166-a97acbaef506.png">
<img width="765" alt="EQN3" src="https://user-images.githubusercontent.com/48260444/209003384-78843ae5-9b65-405e-ad33-07d1d83971b4.png">

1. Acceleration of a particle in the x, y, and z direction due to surrounding masses. This is the only dynamic equation necessary for our FDA approximation.

2. The second-order centred FDA formula which approximates the second derivative of a function.

3. Combining the above two formulas and rearranging the result gives us the desired FDA.

To simplify calculations, we non-dimensionalize by setting the gravitational constant, G, equal to 1. It is also necessary to have two known initial conditions: the initial position and initial velocity of all particles.
