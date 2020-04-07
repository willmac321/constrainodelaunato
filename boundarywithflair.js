import Boundary from './boundary'
import { distLineAndPoint, getEdges, intersect, sortHeap } from './helpers'

export default class BoundaryExtra extends Boundary {
  constructor (arr, k = 3, maxDist = Infinity) {
    super(arr, k, maxDist)
    //    this.cPoints = []
    this.origCoordsLen = arr.length
    this.intersectingLineSegs = []
    this.maxDist = maxDist
  }

  /**
   * addPoints
   * use final k value from concave boundary for point search in order to
   * @param {Array} parentArr Array of coordinate cloud used to interpolate boundary points
   * @param {Object} constructed delaunator object, need this for edges and triangles
   * @param {Integer} dist Max distance to point to trigger interpolation, only one of two points in line segment has to meet this criteria
   */
  addPoints (parentArr, delaunator, dist) {
    this.k = 3
    const edges = getEdges(delaunator)
    // get all intersecting lines to the hull line seg
    for (let p = 0; p < this.hull.length - 1; p++) {
      const h = this.subset([this.hull[p], this.hull[p + 1]])
      const seg = { x0: h[0], y0: h[1], x1: h[2], y1: h[3] }
      const temp = this.getIntersectingLines(seg, edges, parentArr, dist).reverse()
      const ind = sortHeap(temp.map((m) => [m[m.length - 1].x, m[m.length - 1].y]).flat(), [...Array(temp.length).keys()].map((i) => i * 2), 'euclid', [seg.x0, seg.y0])
      this.intersectingLineSegs.push(this.hull[p])
      const c = this.coords.length
      this.coords = this.coords.concat(temp.map((m) => [m[m.length - 1].x, m[m.length - 1].y]).flat())
      this.intersectingLineSegs = this.intersectingLineSegs.concat(ind.map((m) => m + c))
    }
    this.intersectingLineSegs.push(this.intersectingLineSegs[0])

    this.hull = this.intersectingLineSegs
    this.clean(this.hull)
  }

  /**
   * getIntersectingLines
   *
   * @param {Object} lineSeg {x0, y0, x1, y1} Object delineating line segment of boundary or line to find intersecting edges with
   * @param {Array} indexArr Array of indices of edges, these are in pairs where A[n] is index of first point and A[n+1] is index of second point of edge
   * @param {Array} coords Actual coordinate array
   * @param {Float} dist Cutoff distance for intersecting line segments
   *
   * @return {Object} array of an array that contains each point pair index and the x, y coord for intersection format: [index0, index1, {x2, y2}]
   */
  getIntersectingLines (lineSeg, indexArr, coords, dist, opt = null) {
    const pntAndItsArr = []
    // iterate over point pairs
    for (let i = 0; i < indexArr.length; i += 2) {
      // first check if distance from either point to perpendicular of line seg is less than dist
      const point = { x: coords[indexArr[i]], y: coords[indexArr[i] + 1] }
      const d = distLineAndPoint(lineSeg, point)
      if (d <= dist) {
        if (i + 1 > indexArr.length - 1) {
          i = -1
        }
        const coordSeg = { x0: point.x, y0: point.y, x1: coords[indexArr[i + 1]], y1: coords[indexArr[i + 1] + 1] }
        const its = intersect(lineSeg, coordSeg, false)
        if (isFinite(its.x)) {
          // console.log(d, its, lineSeg, coordSeg, i)
          // this.cPoints.push(indexArr[i], indexArr[i + 1])
          if (opt) {
            opt([indexArr[i], indexArr[i + 1], its])
          }
          pntAndItsArr.push([indexArr[i], indexArr[i + 1], its])
        }
      }
    }
    return pntAndItsArr
  }

  get k () {
    return super.k
  }

  set k (v) {
    super.k = v
  }
}
