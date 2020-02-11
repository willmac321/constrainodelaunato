import Delaunator from 'delaunator'
import BoundaryExtra from './boundarywithflair'
import { maximumPointX } from './helpers'

/**
 * ConstrainoDelaunato
 *
 * @class
 * @classdesc ConstrainoDelaunato
 */
export default class ConstrainoDelaunato {
  /**
   * creates a delaunator object for the larger coord point cloud, and any smalle concave boundaries and delaunator objects for holes/boundaries supplied
   * @constructor
   *
   * @param {Array} coords Coordinate cloud, can be 2D or 1D, prefer 1D of type [x0, y0, x1, y1, ... xN, yN]
   * @param {Integer} k lower bound for point selection in k grouping - minimum possible value is 3 - you have to make a polygon
   * @param {Integer} dist - distance for adding points along boundary, distance between line segment perpindicular to either point of triangle segment; used for point interpolation along a boundary
   * @param {Integer} distSelectionLimit - distance to limit selection of candidate points for concave boundary creation - useful if there is a whole on edge of points that is not being acknowledged by algorithm due to uniform point spacing or something like that; used during concave boundary creation
   * @param {Array} ...boundaries Point clouds of holes in coords, stored in array boundary for concave boundaries and boundedDelaunator for created delaunator objects
   */
  constructor (coords, k, dist, distSelectionLimit ...boundaries) {
    // k is the k-nearest neighbor selection
    // if coords are 2D
    if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
      coords = coords.flat()
    } else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
      return
    }

    this._delaunator = new Delaunator(coords)
    this._boundaries = []
    this.boundedDelaunators = []

    for (let boundary of boundaries) {
      if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
        boundary = boundary.flat()
      }
      this._boundaries.push(new BoundaryExtra(boundary, k))
      this._boundaries[this._boundaries.length - 1].addPoints(coords, this._delaunator, dist)
      this.boundedDelaunators.push(this.setTrianglesInsideBound(this._boundaries[this._boundaries.length - 1]))
    }
    this._boundaries.push(new BoundaryExtra(coords, k, distSelectionLimit))
    this._boundaries[this._boundaries.length - 1].addPoints(coords, this._delaunator, dist)
    this.boundedDelaunators.push(this.setTrianglesInsideBound(this._boundaries[this._boundaries.length - 1]))
  }

  /**
   * setTrianglesInsideBound
   *
   * Function used to clip coords to inside of boundary or hole
   *
   * @param {BoundaryExtra} boundary boundary extra object
   */
  setTrianglesInsideBound (boundary) {
    let coords = []
    const outIndex = []
    const index = [...this.delaunator.coords.keys()].filter((i) => i % 2 === 0)
    const maxX = maximumPointX(this.delaunator.coords, index)
    let i = 0
    for (const e of index) {
      const point = { x: this.delaunator.coords[e], y: this.delaunator.coords[e + 1] }
      if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
        outIndex.push(i)
        coords.push(point.x, point.y)
        i += 2
      }
    }

    coords = coords.concat(boundary.subset(boundary.hull))
    const rv = new Delaunator(coords)
    const t = []
    for (let e = 0; e < rv.triangles.length / 3; e++) {
      const edgeIndex = e * 3

      let xCoord = 0
      let yCoord = 0
      for (let r = 0; r < 3; r++) {
        xCoord += rv.coords[2 * rv.triangles[r + edgeIndex]]
        yCoord += rv.coords[2 * rv.triangles[r + edgeIndex] + 1]
      }
      const point = { x: xCoord / 3, y: yCoord / 3 }
      if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
        t.push(rv.triangles[edgeIndex], rv.triangles[edgeIndex + 1], rv.triangles[edgeIndex + 2])
      }
    }
    rv.triangles = new rv.triangles.constructor(t)

    return rv
  }

  /**
   * update
   *
   * @param {Array} point x and y coord of point to add the delaunator object
   */
  update (point) {
    const c = this.coords
    for (const p of point.flat()) {
      c.push(p)
    }
    this._delaunator = new Delaunator(c)
  }

  /**
   * coords2D
   *
   * @returns {Array} 2D coordinate array
   */
  get coords2D () {
    const c2D = []
    const c1D = this.coords
    for (let i = 0; i < c1D.length; i += 2) {
      c2D.push([c1D[i], c1D[i + 1]])
    }
    return c2D
  }

  /**
   * coords
   *
   * @returns {Array} 1D coordinate array
   */
  get coords () {
    return this.delaunator.coords
  }

  /**
   * triangles
   *
   * @returns {Array} Index array of delaunator triangles
   */
  get triangles () {
    return this.delaunator.triangles
  }

  /**
   * hull
   *
   * @returns {Array} Array of hull indices
   */
  get hull () {
    return this.delaunator.hull
  }

  /**
   * delaunator
   *
   * @returns {Object} Return delaunator object for parent points
   */
  get delaunator () {
    return this._delaunator
  }

  /**
   * boundaries
   *
   * @returns {Array} concave boundary index array of parent points and hole points includes parent points as final array item
   */
  get boundaries () {
    return this._boundaries
  }

  /**
   * holes
   *
   * @returns {Array} Delaunator object array for all hole/boundary points supplied, includes parent points as final array item
   */
  get holes () {
    return this.boundedDelaunators
  }
}
