import { euclid, intersect, maximumPointX, maximumPointY, minimumPointY, minimumPointX, slope, swap, sortHeap } from './helpers'

var counter = 0

export default class Boundary {
  constructor (arr, k = 3) {
    this.k = k
    this.coords = arr.slice()
    this.index = [...this.coords.keys()].filter((i) => i % 2 === 0)
    this.index = this.clean(this.index)
    this.center = this.calcCenter()
    this.minY = minimumPointY(this.coords, this.index)
    this.maxY = maximumPointY(this.coords, this.index)
    this.minX = minimumPointX(this.coords, this.index)
    this.maxX = maximumPointX(this.coords, this.index)
    this.maxD = Math.sqrt(Math.pow(this.maxX.x - this.minX.x, 2) + Math.pow(this.maxY.y - this.minY.y, 2))
    this.maxR = this.maxR / 2
    this.offsetAngle = 1

    this.cPoints = []

    this.ray = null
    this.hull = this.findConcaveHull(k)
  }

  /**
   * findConcaveHull
   *
   * @param {Integer} k Starting point cloud count for points being used for selection
   * @returns {Array} Array of indices of X values for the sorted concave hull
   */
  findConcaveHull (k) {
    // alt index is sorted to minX value
    const index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y])
    const hull = this.concave(index, k)
    // hull = this.sortHeapAndClean(this.coords, hull, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y])
    // hull.push(hull[0])
    return hull
  }

  concave (index, k) {
    // k nearest neighbor babbbbyyyy
    // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
    // double check arr is sorted and clean
    // also sort it so all points are in order from some min point  on the xy plane
    const stopVal = Infinity // 200 // 86 // Infinity // 76 // Infinity // and beyond
    const oldIndex = index.slice()
    // console.log('new k', k)
    if (index.length < 3) {
      console.log('len less than 3')
      return null
    } else if (k > index.length - 1) {
      console.log(counter)
      console.log('k is too big')
      return null
    } else if (index.length === 3) {
      console.log('len 3')
      return index
    }

    let kk = Math.min(Math.max(k, 3), index.length - 1)
    // i is a pointer to the relative index not a loc in this.coords
    // so, index of that index gives a this.coords pointer
    const firstPointIndex = minimumPointY(this.coords, index).i
    const firstPoint = { i: firstPointIndex, coord: index[firstPointIndex] }
    let currentPoint = firstPoint.coord
    const hull = [firstPoint.coord]
    // why is step init to 2?
    // Because the paper was written in Matlab....
    let step = 1
    // each index value can only be used once so this is ok
    index.splice(firstPoint.i, 1)
    while ((currentPoint !== firstPoint.coord || step === 1) && (index.length > 0)) {
      counter++
      if (step === 4) {
        index.push(firstPoint.coord)
      }
      // find nearest neighbors
      const kNearestPoints = this.nearestPoints(index, currentPoint, kk)
      // descending order 'right-hand' turn x and y min are top left on js canvas in webpage
      const cPoints = this.sortByAngle(kNearestPoints, currentPoint, hull[hull.length - 2])
      // if (cPoints.indexOf(firstPoint.coord) > -1) {
      //   console.log(cPoints)
      // }
      let its = true
      let i = -1
      while (its && i < cPoints.length - 1) {
        // This is so that when the first point is added to the end of the hull, it doesn't get used to check for intersections
        let lastPoint = 0
        if (cPoints[i] === firstPoint.coord) {
          // console.log('back to first', firstPoint)
          lastPoint = 1
        }
        let j = 1
        its = false
        while (!its && j < hull.length - lastPoint) {
          const l = {
            x0: this.coords[hull[step - 1]],
            y0: this.coords[hull[step - 1] + 1],
            x1: this.coords[cPoints[i + 1]],
            y1: this.coords[cPoints[i + 1] + 1]
          }
          const p = {
            x0: this.coords[hull[step - j]],
            y0: this.coords[hull[step - j] + 1],
            x1: this.coords[hull[step - 1 - j]],
            y1: this.coords[hull[step - 1 - j] + 1]
          }
          // the endpoint of one line segment is always intersecting the endpoint of a connected line segment, how to ignore this intersection?
          const ints = intersect(p, l, true)
          const endpointsMatch = (p.x0 === l.x0 && p.y0 === l.y0)
          const isClose = (cPoints[i + 1] === firstPoint.coord) && (p.x1 === l.x1 && p.y1 === l.y1)
          // (p.x0 !== l.x0 && p.y0 !== l.y0) ||
          // if (l.x0 === 221 && l.y0 === 90) {
          //   console.log(l, p, ints, isFinite(ints.x), (p.x1 === l.x0 && p.y1 === l.y0))
          // }
          if (isFinite(ints.x) && !endpointsMatch && !isClose) {
            its = true
          }
          j++
        }
        i++
      }
      this.cPoints = cPoints.slice()
      // this.cPoints.splice(i, 1)

      if (its) {
        return this.concave(oldIndex, ++kk)
      }
      currentPoint = cPoints[i]
      hull.push(currentPoint)

      if (counter > stopVal) {
        return hull // .concat(cPoints)
      }
      index.splice(index.indexOf(currentPoint), 1)
      step++
    }
    let allInside = true
    for (const i of index) {
      allInside = this.pointInOrOut(
        [this.coords[i], this.coords[i + 1]],
        hull, 0)
      if (!allInside) {
        break
      }
    }
    if (!allInside) {
      return this.concave(oldIndex, ++kk)
    }
    this.k = kk
    return hull
  }

  sortByAngle (kNearestPoints, currentPoint, lastPoint) {
    if (!lastPoint || lastPoint === currentPoint) {
      lastPoint = [this.maxX.x + 10, this.coords[currentPoint + 1]]
    } else {
      lastPoint = [this.coords[lastPoint], this.coords[lastPoint + 1]]
    }
    this.ray = {
      x0: lastPoint[0],
      y0: lastPoint[1],
      x1: this.coords[currentPoint],
      y1: this.coords[currentPoint + 1]
    }
    const currentPointArr = [this.coords[currentPoint], this.coords[currentPoint + 1]]
    // cant use max or min value for first point, the reference point needs to be the last point in the hull in order to get the angle sorting right
    const rv = sortHeap(this.coords, kNearestPoints, 'polar', lastPoint, currentPointArr).slice()
    // if two points are on the same line eq as current point, currently the further one is considered a 'closer angle', perform swap of these coords below

    let lastSlope
    let lastDist
    // if two points relative to each other are in line
    // Issue here when 3 points line up and one is segment from origin
    for (let k = 0; k < rv.length; k++) {
      let lastPoint = [this.coords[rv[k - 1]], this.coords[rv[k - 1] + 1]]
      if (k === 0) {
        lastPoint = currentPointArr
      }
      const newPoint = [this.coords[rv[k]], this.coords[rv[k] + 1]]
      const newSlope = slope(lastPoint, newPoint)
      const newDist = euclid(currentPointArr, newPoint)

      if (!isNaN(lastSlope) && !isNaN(lastDist) && (Math.abs(newSlope) === Math.abs(lastSlope) || (newSlope === Infinity && lastSlope === -Infinity)) && newDist < lastDist) {
        // flipflop the two points in array order if the slopes are the same
        // sort by euclid instead of straight swap
        swap(rv, k, k - 1)
        lastDist = euclid(currentPointArr, [this.coords[rv[k]], this.coords[rv[k] + 1]])
      } else {
        lastDist = newDist
      }

      lastSlope = slope(lastPoint, newPoint)
    }

    return rv
  }

  nearestPoints (index, cP, kk) {
    const currentPoint = [this.coords[cP], this.coords[cP + 1]]
    index = sortHeap(this.coords.slice(), index.slice(), 'euclid', currentPoint)
    const rv = []
    let lastSlope
    let lastDist
    kk = Math.min(kk, index.length - 1)
    let i = 0
    let c = 0
    while (c < kk) {
      const newSlope = slope(currentPoint, [this.coords[index[i]], this.coords[index[i] + 1]])
      const newDist = euclid(currentPoint, [this.coords[index[i]], this.coords[index[i] + 1]])

      if (newDist === lastDist) {
        rv.push(index[i])
      } else if (c === 0 || newSlope !== lastSlope) {
      // } else if ( !lastSlope || newSlope !== lastSlope) {
        rv.push(index[i])
        c++
      }
      i++
      if (i > index.length - 1) {
        return rv
      }
      lastSlope = newSlope
      lastDist = newDist
    }
    return rv
  }

  sortHeapAndClean (arr, ind, criteria, minPoint, centerPoint) {
    ind = sortHeap(arr.slice(), ind.slice(), criteria, minPoint, centerPoint)
    ind = this.clean(ind)
    return ind
  }

  clean (index) {
    // there has to be a better way to do this
    // On^2  urrrgh
    const itRem = index.length

    let count = 0
    const duplicates = []
    for (const item of index) {
      for (let i = 0; i < index.length; i++) {
        if (this.coords[index[i]] === this.coords[item] &&
          this.coords[index[i] + 1] === this.coords[item + 1] &&
          count !== i) {
          let pass = true
          for (const t of duplicates) {
            if (this.coords[index[i]] === this.coords[t[0]] && this.coords[index[i] + 1] === this.coords[t[0] + 1]) {
              pass = false
              t[1]++
            }
          }
          if (pass) { duplicates.push([item, 0]) }
          break
        }
      }
      count++
    }
    const newIndex = []
    for (let i = 0; i < index.length; i++) {
      let pass = true
      for (const item of duplicates) {
        if (this.coords[index[i]] === this.coords[item[0]] &&
          this.coords[index[i] + 1] === this.coords[item[0] + 1] &&
          item[1] > 0) {
          item[1]--
          pass = false
        }
      }
      if (pass) { newIndex.push(index[i]) }
    }
    // console.log('items removed: ' + (itRem - newIndex.length))
    return newIndex
  }

  calcCenter () {
    const p = { x: 0, y: 0 }

    for (let i = 0; i < this.coords.length; i += 2) {
      p.x += this.coords[i]
      p.y += this.coords[i + 1]
    }
    p.x /= (this.coords.length / 2)
    p.y /= (this.coords.length / 2)
    return p
  }

  /**
   * pointInOrOut
   *
   * checks if a point is in or out of the polygon, if a vertex is encountered at an angle, will rotate the ray around the point until a weird intersection is not found -- weird here is the intersection of two linesegments of the boundary
   *        - pretty much avoid address intersect boundary line segments by rotating the array
   *
   * @param {Array} point point to use as center point for ray to test for simple closed polygon
   * @param {Array} index index array of the boundary points
   * @param {Integer} angle Incremental angle to change the ray by if a weird vertex is encountered, will quit at 360 degree rotation, hopefully theres an angle there, if not will just return a value
   * @returns {Boolean} True if the point is inside, false if the point is outside the polygon
   */
  pointInOrOut (point, index, angle) {
    // assume ray going to + infinity on x plane here just making assumption that it extends 1000 units past whatever the minimum x value is in the boundary
    const x1 = point[0] + this.maxD * Math.cos(angle * Math.PI / 180)
    const y1 = point[1] + this.maxD * Math.sin(angle * Math.PI / 180)
    const p = {
      x0: point[0], y0: point[1], x1: x1, y1: y1
    }
    this.ray = p
    // lets use non-zero winding number rule
    let windingNum = 0
    let last = { x: Infinity, y: Infinity }

    for (let i = 0; i < index.length; i++) {
      const l = {
        x0: this.coords[index[i]],
        y0: this.coords[index[i] + 1],
        x1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1]],
        y1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1] + 1]
      }
      const inters = intersect(p, l, true)
      if (isFinite(inters.x)) {
        // fails on corner intersection sometimes, so if a corner is intersect, rotate the ray being tested by 1 deg and start over
        const testCond = Math.round(inters.x * 1000000) === last.x && Math.round(inters.y * 1000000) === last.y
        if (l.y1 - l.y0 > 0 && !testCond) {
          windingNum++
        } else if (l.y1 - l.y0 < 0 && !testCond) {
          windingNum--
        } else if (testCond) {
          if (angle < 360) {
            return this.pointInOrOut(point, index, ++angle)
          }
        }
        last = { x: Math.round(inters.x * 1000000), y: Math.round(inters.y * 1000000) }
      }
    }
    return Math.abs(windingNum) !== 0
  }

  printPoints (xIndex) {
    const p = []
    for (const i of xIndex) {
      p.push(this.coords[i], this.coords[i + 1])
    }
    console.log(p)
  }

  get hullCoords () {
    return this.subset(this.hull)
  }

  subset (indices) {
    const rv = []
    for (const i of indices) {
      rv.push(this.coords[i], this.coords[i + 1])
    }
    return rv
  }

  get coords2D () {
    const newArr = []
    const arr = this.sortedCoords
    while (arr.length) newArr.push(arr.splice(0, 2))
    return newArr
  }

  get sortedCoords () {
    const newArr = []
    for (const i of this.hull) {
      newArr.push(this.coords[i], this.coords[i + 1])
    }
    return newArr
  }
}
