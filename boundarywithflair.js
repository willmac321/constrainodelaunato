import Boundary from './boundary'
import { distLineAndPoint, intersect } from './helpers'

export default class BoundaryExtra extends Boundary {
  constructor (arr, k = 3) {
    super(arr, k)
    this.cPoints = []
  }

  /**
   * addPoints
   * use final k value from concave boundary for point search in order to
   * @param {Array} parentArr Array of coordinate cloud used to interpolate boundary points
   * @param {Array} parentIndex Sorted index array of coordinate cloud used to interpolate boundary points
   * @param {Integer} dist Max distance to point to trigger interpolation, only one of two points in line segment has to meet this criteria
   */
  addPoints (parentArr, parentIndex, dist) {
    this.k = 3
    // get all intersecting lines to the hull line seg
    for (let p = 0; p < this.hull.length - 1; p++) {
      const h = this.subset([this.hull[p], this.hull[p + 1]])
      const seg = { x0: h[0], y0: h[1], x1: h[2], y1: h[3] }
      this.getIntersectingLines(seg, parentIndex, parentArr, dist)
    }
  }

  getIntersectingLines (lineSeg, indexArr, coords, dist) {
    for (let [i, t] of indexArr.entries()) {
      // first check if distance from either point to perpendicular of line seg is less than dist
      const point = { x: coords[t], y: coords[t + 1] }
      const d = distLineAndPoint(lineSeg, point)
      if (d <= dist) {
        if (i + 1 > indexArr.length - 1) {
          i = -1
        }
        const coordSeg = { x0: point.x, y0: point.y, x1: coords[indexArr[i + 1]], y1: coords[indexArr[i + 1] + 1] }
        console.log(coordSeg)
        return
        const its = intersect(lineSeg, coordSeg)
        if (isFinite(its.x)) {
          console.log(d, its)
        }
      }
    }
  }

  interpolate (parentArr) {

  }

  get k () {
    return super.k
  }

  set k (v) {
    super.k = v
  }
}
