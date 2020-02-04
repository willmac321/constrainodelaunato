import Delaunator from 'delaunator'
import Boundary from './boundarywithflair'
import { maximumPointX } from './helpers'

export default class ConstrainoDelaunato {
  constructor (coords, k, ...boundaries) {
    // k is the k-nearest neighbor selection
    // if coords are 2D
    if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
      coords = coords.flat()
    } else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
      return
    }

    this.delaunator = new Delaunator(coords)
    this.boundaries = []
    this.boundedDelaunators = []

    for (let boundary of boundaries) {
      if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
        boundary = boundary.flat()
      }
      if (boundary) {
        this.boundaries.push(new Boundary(boundary, k))
      } else {
        this.boundaries.push(new Boundary(coords, k))
      }
      this.boundaries[this.boundaries.length - 1].addPoints(coords, this.delaunator, 10)
      this.boundedDelaunators.push(this.setTrianglesInsideBound(this.boundaries[this.boundaries.length - 1]))
    }

    this.boundary = this.boundaries[this.boundaries.length - 1]
    this.boundedDelaunator = this.boundedDelaunators[this.boundedDelaunators.length - 1]
  }

  setTrianglesInsideBound (boundary) {
    let coords = []
    const outIndex = []
    const index = [...this.delaunator.coords.keys()].filter((i) => i % 2 === 0)
    const maxX = maximumPointX(this.delaunator.coords, index)
    let i = 0
    for (const e of index) {
      const point = { x: this.delaunator.coords[e], y: this.delaunator.coords[e + 1] }
      // if (point.x === 59 && point.y === 80) {
      //   console.log(boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10))
      // }
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
