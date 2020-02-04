import Delaunator from 'delaunator'
import Boundary from './boundarywithflair'
import { getEdges, maximumPointX } from './helpers'

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
    } else {
      this.boundary = new Boundary(coords, k)
    }

    this.delaunator = new Delaunator(coords)

    this.boundary.addPoints(coords, this.delaunator, 10)
    this.boundedDelaunator = this.setTrianglesInsideBound(this.boundary)
    this.delaunator = this.boundedDelaunator
  }

  setTrianglesInsideBound (boundary) {
    let coords = []
    const outIndex = []
    const index = [...this.delaunator.coords.keys()].filter((i) => i % 2 === 0)
    const maxX = maximumPointX(this.delaunator.coords, index)
    let i = 0
    for (const e of index) {
      const point = { x: this.delaunator.coords[e], y: this.delaunator.coords[e + 1] }
      if (point.x === 59 && point.y === 80) {
        console.log(boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10))
      }
      if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
        outIndex.push(i)
        coords.push(point.x, point.y)
        i += 2
      }
    }

    coords = coords.concat(boundary.subset(boundary.hull))
    const rv = new Delaunator(coords)
    // TODO have to remove triangs from the convex bound created
    return rv
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
