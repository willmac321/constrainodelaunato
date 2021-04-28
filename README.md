# Constraino Delaunato
### your a wizard harry
<img src="https://gitlab.com/illmac321/constrainodelaunato/-/raw/master/static/wizard.jpg" alt="Tada" width="300"/>

## About

<div style="align:center">
<img src="https://gitlab.com/willmac321/constrainodelaunato/-/raw/master/static/d4.png" alt="Constraino Delaunato Example" width = "300"/>
</div>

Based on/using delaunator as basis for delaunay triangulations, this module uses the algorithm proposed in this paper, 
[[CONCAVE HULL: A K-NEAREST NEIGHBOURS APPROACH FOR THE COMPUTATION OF THE REGION OCCUPIED BY A SET OF POINTS](https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf)],
in order to create a concave hull around a point cloud.

- Additionally, the module takes in point sets as additional optional arguments.  When these are supplied, the a seperate concave boundary is created for the boundary points.  
    - These points also have their own triagulation created by taking available points from the source data and clipping that area out to match the created concave boundary.

## Example 
- A delaunay triangulation with a new set of points overlayed in black - the black line is to help with visibility of the points, the points are at vertices of the polygon:
<div style="align:center"><img src="https://gitlab.com/willmac321/constrainodelaunato/-/raw/master/static/d1.png" alt="Constraino Delaunato Ex 1" width = "400"/></div>

- The same set of points but with a concave boundary created around it in blue, the arrows indicate a ccw direction.  (CCW from the origin (0,0) which is at the top left point of this image)  The concave algorithm follows the left hand rule to determine the most suitable next point in the k grouping.  
    - At this point, additional points are added along the boundary for intersections with the main triangulation.  This follows a minimum distance logic, where the points are not considered for an intersection with the boundary unless one of the points perpindicular line intersects the boundary and is less than or equal to that minimum distance from that intersection. [[Further Details](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line)]
<div style="align:center">
<img src="https://gitlab.com/willmac321/constrainodelaunato/-/raw/master/static/d2.png" alt="Constraino Delaunato Ex 2" width = "400"/>
</div>

- Finally the concave boundary of the inner points are used to clip the original points into a seperate delaunay object.  This object then has any triangle outside of the new concave boundary removed.  It is done in this order to in order to minimize the size of the initial delaunay object and also to include the boudnary intersection points in the triangulation calculation.  The new triangulation is shown overlayed in blue.
<div style="align:center">
<img src="https://gitlab.com/willmac321/constrainodelaunato/-/raw/master/static/d3.png" alt="Constraino Delaunato Ex 3" width = "400"/>
</div>

## Run

### For Dev - if cloned from git:
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
<a name="ConstrainoDelaunato"></a>

## ConstrainoDelaunato
ConstrainoDelaunato

**Kind**: global class  

* [ConstrainoDelaunato](#ConstrainoDelaunato)
    * [new ConstrainoDelaunato(coords, k, dist)](#new_ConstrainoDelaunato_new)
    * [.coords2D](#ConstrainoDelaunato+coords2D) ⇒ <code>Array</code>
    * [.coords](#ConstrainoDelaunato+coords) ⇒ <code>Array</code>
    * [.triangles](#ConstrainoDelaunato+triangles) ⇒ <code>Array</code>
    * [.hull](#ConstrainoDelaunato+hull) ⇒ <code>Array</code>
    * [.delaunator](#ConstrainoDelaunato+delaunator) ⇒ <code>Object</code>
    * [.boundaries](#ConstrainoDelaunato+boundaries) ⇒ <code>Array</code>
    * [.holes](#ConstrainoDelaunato+holes) ⇒ <code>Array</code>
    * [.setTrianglesInsideBound(boundary)](#ConstrainoDelaunato+setTrianglesInsideBound)
    * [.update(point)](#ConstrainoDelaunato+update)

<a name="new_ConstrainoDelaunato_new"></a>

### new ConstrainoDelaunato(coords, k, dist)
creates a delaunator object for the larger coord point cloud, and any smalle concave boundaries and delaunator objects for holes/boundaries supplied


| Param | Type | Description |
| --- | --- | --- |
| coords | <code>Array</code> | Coordinate cloud, can be 2D or 1D, prefer 1D of type [x0, y0, x1, y1, ... xN, yN] |
| k | <code>Integer</code> | lower bound for point selection in k grouping - minimum possible value is 3 - you have to make a polygon |
| dist | <code>Integer</code> | distance for adding points along boundary, distance between line segment perpindicular to either point of triangle segment |
| ...boundaries | <code>Array</code> | Point clouds of holes in coords, stored in array boundary for concave boundaries and boundedDelaunator for created delaunator objects |

<a name="ConstrainoDelaunato+coords2D"></a>

### constrainoDelaunato.coords2D ⇒ <code>Array</code>
coords2D

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - 2D coordinate array  
<a name="ConstrainoDelaunato+coords"></a>

### constrainoDelaunato.coords ⇒ <code>Array</code>
coords

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - 1D coordinate array  
<a name="ConstrainoDelaunato+triangles"></a>

### constrainoDelaunato.triangles ⇒ <code>Array</code>
triangles

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - Index array of delaunator triangles  
<a name="ConstrainoDelaunato+hull"></a>

### constrainoDelaunato.hull ⇒ <code>Array</code>
hull

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - Array of hull indices  
<a name="ConstrainoDelaunato+delaunator"></a>

### constrainoDelaunato.delaunator ⇒ <code>Object</code>
delaunator

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Object</code> - Return delaunator object for parent points  
<a name="ConstrainoDelaunato+boundaries"></a>

### constrainoDelaunato.boundaries ⇒ <code>Array</code>
boundaries

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - concave boundary index array of parent points and hole points includes parent points as final array item  
<a name="ConstrainoDelaunato+holes"></a>

### constrainoDelaunato.holes ⇒ <code>Array</code>
holes

**Kind**: instance property of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  
**Returns**: <code>Array</code> - Delaunator object array for all hole/boundary points supplied, includes parent points as final array item  
<a name="ConstrainoDelaunato+setTrianglesInsideBound"></a>

### constrainoDelaunato.setTrianglesInsideBound(boundary)
setTrianglesInsideBound

Function used to clip coords to inside of boundary or hole

**Kind**: instance method of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  

| Param | Type | Description |
| --- | --- | --- |
| boundary | <code>BoundaryExtra</code> | boundary extra object |

<a name="ConstrainoDelaunato+update"></a>

### constrainoDelaunato.update(point)
update

**Kind**: instance method of [<code>ConstrainoDelaunato</code>](#ConstrainoDelaunato)  

| Param | Type | Description |
| --- | --- | --- |
| point | <code>Array</code> | x and y coord of point to add the delaunator object |


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
