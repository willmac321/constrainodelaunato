import Delaunator from 'delaunator'
import Boundary from './boundarywithflair'

export default class ConstrainoDelaunato {
  constructor (coords, boundary, k) {
    // k is the k-nearest neighbor selection
    // if coords are 2D
    if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
      coords = coords.flat()
    } else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
      return
    }
    if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
      boundary = boundary.flat()
    }
    if (boundary) {
      this.boundary = new Boundary(boundary, k)
      //coords = coords.concat(this.boundary.hullCoords)
    }
    this.delaunator = new Delaunator(coords)
    this.boundary.addPoints(coords, this.delaunator, 10)
    // this.pointInOrOut([1,1]);
  }

  update (point) {
    const c = this.coords
    for (const p of point.flat()) {
      c.push(p)
    }
    this.delaunator = new Delaunator(c)
  }

  get coords2D () {
    const c2D = []
    const c1D = this.coords
    for (let i = 0; i < c1D.length; i += 2) {
      c2D.push([c1D[i], c1D[i + 1]])
    }
    return c2D
  }

  get coords () {
    return this.delaunator.coords
  }

  get triangles () {
    return this.delaunator.triangles
  }

  get concaveHullCoords () {
    return this.boundary.hullCoords
  }

  get hull () {
    return this.delaunator.hull
  }

  get bound () {
    return this.boundary
  }
}
