# Constraino Delaunato
### your a wizard harry
![tada](https://i.redd.it/qc8idcyojj231.jpg)

## About
Based on/using delaunator as basis for delaunay triangulations, this module uses the algorithm proposed in this paper, 
[![CONCAVE HULL: A K-NEAREST NEIGHBOURS APPROACH FOR THE COMPUTATION OF THE REGION OCCUPIED BY A SET OF POINTS](https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf)],
in order to create a concave hull around a point cloud.

- Additionally, the module takes in othe point sets as additional optional arguments.  When these are supplied, the a seperate concave ooundary is created for the boundary points.  
    - These points also have their own triagulation created by taking available points from the source data and clipping that area out to match the created concave boundary.

## Example 
- A delaunay triangulation with a new set of points overlayed in black:
<div style="text-align:center"><img src="/assets/d1.png" alt="Constraino Delaunato Ex 1" width = "400"/></div>

- The same set of points but with a concave boundary created around it in blue, the arrows indicate a ccw direction.  (CCW from the origin (0,0) which is at the top left point of this image)  The concave algorithm follows the left hand rule to determine the most suitable next point in the k grouping.  
    - At this point, additional points are added along the boundary for intersections with the main triangulation.  This follows a minimum distance logic, where the points are not considered for an intersection with the boundary unless one of the points perpindicular line intersects the boundary and is less than or equal to that minimum distance from that intersection. [![Further Details](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line)]
<div style="text-align:center">
<img src="/assets/d2.png" alt="Constraino Delaunato Ex 2" width = "400"/>
</div>

- Finally the concave boundary of the inner points are used to clip the original points into a seperate delaunay object.  This object then has any triangle outside of the new concave boundary removed.  It is done in this order to in order to minimize the size of the initial delaunay object and also to include the boudnary intersection points in the triangulation calculation.  The new triangulation is shown overlayed in blue.
<div style="text-align:center">
<img src="/assets/d3.png" alt="Constraino Delaunato Ex 3" width = "400"/>
</div>

## Run

### For Dev:
- uses rollup for minimizing
```bash
npm install [package name]
npm start
```

### For not dev
```js
//ES import
import ConstrainoDelaunato as CD from 'constrainodelaunato'

//this is it, it will do the hole triangle calc on its own
let rictusempra = new CD(points, 3, holeBoundary)

// can add multiple holes
let tarantallegra = new CD(points, 3, hole1, hole2, holeN)

// hole boundaries and hole delaunay objects are stored as seperate arrays
console.log(tarantallegra.boundaries[0])
console.log(tarantallegra.boundedDelaunators[0])
```
- The bounded delaunator objects for the holes are useable as delaunay objects

### API

#### Everything below is done automatically - but just in case...

#### BoundaryExtra

```js
// 3 is k starting value
let b = new BoundaryExtra(points, 3)
b.addPoints(parentCoords, parentDelaunatorObject, distance)
```

#### below if for boundary object - should use boundaryExtra instead though
```js
// for Object.boundaries 3 is k starting value
let b = new Boundary(points, 3)

//below is already done on object creation
const hull = b.findConcaveHull(k)// -> returns: sorted array of indices of x values for coord array

//subset returns the subset of coordinates from an index array
console.log(b.subset[hull])    

//other accessors
b.coords2D
b.hullCoords
b.printPoints
```
