import Delaunator from 'delaunator'
import { heapSort, intersect, dotProduct } from './helpers'

// const test = [10, 6, 3, 4, 7, 1, 2, 5]

var counter = 0

class Boundary {
  constructor (arr) {
    this.coords = arr.slice()
    this.index = [...this.coords.keys()].filter((i) => i % 2 === 0)
    this.index = this.clean(this.index)
    this.center = this.calcCenter()
    this.minY = minimumPointY(this.coords, this.index)
    this.minX = minimumPointX(this.coords, this.index)
    this.maxX = maximumPointX(this.coords, this.index)

    this.ray = null
    this.hull = []

    // center of test from html is not inside boundary
    // this point is though
    this.testFunctions()
  }

  testFunctions () {
    this.pointInOrOut([this.center.x, this.center.y], this.index, this.minX.x - 10)
    this.pointInOrOut([this.minX.x + 1000, this.minX.y], this.index, this.minX.x - 10)
    // console.log(this.pointInOrOut([180, 100], this.index))
    this.index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minY.x, this.minY.y], [this.minX.x, this.minY.y])
    this.findHull(3)
  }

  findHull (k) {
    // alt index is sorted to minX value
    this.index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minX.x, this.minX.y], [this.minY.x, this.minY.y])
    // this.index = index
    this.index = this.concave(this.index.slice(), k)
    // this.index = this.index.filter((i) => i !== undefined)
    console.log(this.index)
    // this.sortHeapAndClean(this.coords, this.index, 'polar', [this.center.x, this.center.y])
  }

  concave (index, k) {
    // k nearest neighbor babbbbyyyy
    // https://towardsdatascience.com/the-concave-hull-c649795c0f0f
    // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
    // double check arr is sorted and clean
    // also sort it so all points are in order from some min point  on the xy plane
    const stopVal = Infinity
    const oldIndex = index.slice()
    console.log('new k', k)
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

    // TODO why is step init to 2?
    // Because the paper was written in Matlab....
    let step = 1

    // each index value can only be used once so this is ok
    // console.log(firstPoint)
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
      // console.log('cPoints', cPoints.slice(), kNearestPoints, this.subset(cPoints))

      let its = true
      let i = -1
      // console.log('here', cPoints[i + 1], firstPoint)

      while (its && i < cPoints.length - 1) {
        // TODO quit early to check drawing

        // This is so that when the first point is added to the end of the hull, it doesn't get used to check for intersections
        let lastPoint = 0
        if (cPoints[i] === firstPoint.coord) {
          lastPoint = 1
        }

        // console.log(`last point ${lastPoint}`)

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
            x0: this.coords[hull[step - 1 - j]],
            y0: this.coords[hull[step - 1 - j] + 1],
            x1: this.coords[hull[step - j]],
            y1: this.coords[hull[step - j] + 1]
          }
          // console.log(l, p)
          // console.log('hull', this.subset(hull))
          its = isFinite(intersect(p, l).x)
          j++
        }
        i++
      }
      // console.log(its)
      if (its) {
        console.log('intersection found at k ', k, its)
        return this.concave(oldIndex, ++kk)
      }
      currentPoint = cPoints[i]
      hull.push(currentPoint)
      if (counter > stopVal) {
        return hull
      }
      // previousAngle = dotProduct(
      //   [this.coords[currentPoint], this.coords[currentPoint + 1]],
      //   [this.coords[hull[hull.length - 2]], this.coords[hull[hull.length - 2] + 1]])
      index.splice(index.indexOf(currentPoint), 1)
      step++
    }
    let allInside = true
    for (const i of index) {
      console.log(this.coords[i], this.coords[i + 1])
      allInside = this.pointInOrOut(
        [this.coords[i], this.coords[i + 1]],
        hull, this.maxX.x + 10)
      if (!allInside) {
        break
      }
    }
    if (!allInside) {
      console.log('Another time round')
      return hull
      return this.concave(oldIndex, ++kk)
    }
    console.log('made it out')
    return hull
  }

  sortByAngle (kNearestPoints, currentPoint, lastPoint) {
    // TODO does this work as expected?

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

    // cant use max or min value for first point, the reference point needs to be the last point in the hull in order to get the angle sorting right
    return sortHeap(this.coords, kNearestPoints, 'polar', lastPoint, [this.coords[currentPoint], this.coords[currentPoint + 1]])
  }

  nearestPoints (index, cP, kk) {
    // console.log(cP)
    // console.log([this.coords[cP], this.coords[cP + 1]])
    index = sortHeap(this.coords.slice(), index.slice(), 'dist', [this.coords[cP], this.coords[cP + 1]])
    // console.log(index)
    const rv = []
    kk = Math.min(kk, index.length - 1)
    for (let i = 0; i < kk; i++) {
      rv.push(index[i])
    }
    return rv
  }

  sortHeapAndClean (arr, ind, criteria, minPoint, centerPoint) {
    // console.log(this.index, this.coords2D)
    console.log(minPoint, centerPoint)
    ind = sortHeap(arr.slice(), ind.slice(), criteria, minPoint, centerPoint)
    console.log('heap clean res\n', arr, ind)
    ind = this.clean(ind)
    return ind
  }

  clean (index) {
    // TODO there has to be a better way to do this
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
    console.log('items removed: ' + (itRem - newIndex.length))
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

  pointInOrOut (point, index, dir) {
    // assume ray going to + infinity on x plane here just making assumption that it extends 1000 units past whatever the minimum x value is in the boundary
    const p = {
      x0: point[0], y0: point[1], x1: dir, y1: point[1]
    }
    this.ray = p
    // console.log(this.ray)
    // lets use non-zero winding number rule
    let windingNum = 0

    for (let i = 0; i < index.length; i++) {
      const l = {
        x0: this.coords[index[i]],
        y0: this.coords[index[i] + 1],
        x1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1]],
        y1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1] + 1]
      }
      const inters = intersect(p, l, true)
      if (isFinite(inters.x)) {
        if (l.y1 - l.y0 > 0) {
          windingNum++
        } else if (l.y1 - l.y0 < 0) {
          windingNum--
        }
      }
    }
    console.log(windingNum !== 0)
    return Math.abs(windingNum) !== 0
  }

  printPoints (xIndex) {
    const p = []
    for (const i of xIndex) {
      p.push(this.coords[i], this.coords[i + 1])
    }
    console.log(p)
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
    for (const i of this.index) {
      newArr.push(this.coords[i], this.coords[i + 1])
    }
    return newArr
  }
}

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
      coords = coords.concat(this.boundary.coords)
    }
    this.delaunator = new Delaunator(coords)
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

  get hull () {
    return this.delaunator.hull
  }

  get bound () {
    return this.boundary.sortedCoords.flat()
  }
}

function sortHeap (arr, index, criteria, minPoint, centerPoint) {
  // convert point arr to 2d -> easier for me to get my head around sorting

  // minPoint = { x: minPoint.x, y: minPoint.y }
  // builtInSort([minX, minY], newArr);

  heapSort(minPoint, index, arr, index.length, criteria, centerPoint)

  return index
}

function maximumPointX (newArr, index) {
  let ind = 0
  let minY = -Infinity
  let minX = -Infinity
  if (index) {
    for (const [k, p] of index.entries()) {
      if (newArr[p] > minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      } else if (newArr[p + 1] >= minY && newArr[p] >= minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      }
    }
  } else {
    for (let p = 0; p > newArr.length; p++) {
      if (newArr[p] > minX) {
        minX = newArr[p]
        ind = p
      }
    }
  }
  return { x: minX, y: minY, i: ind }
}

function minimumPointY (newArr, index) {
  let ind = 0
  let minY = Infinity
  let minX = Infinity
  if (index) {
    for (const [k, p] of index.entries()) {
      if (newArr[p + 1] < minY) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      }
    }
  } else {
    for (let p = 0; p < newArr.length; p++) {
      if (newArr[p] < minX) {
        minX = newArr[p]
        ind = p
      }
    }
  }
  console.log({ x: minX, y: minY, i: ind })
  return { x: minX, y: minY, i: ind }
}

function minimumPointX (newArr, index) {
  let ind = 0
  let minY = Infinity
  let minX = Infinity
  if (index) {
    for (const [k, p] of index.entries()) {
      if (newArr[p] < minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = k
      }
    }
  } else {
    for (let p = 0; p < newArr.length; p++) {
      if (newArr[p] < minX) {
        minX = newArr[p]
        ind = p
      }
    }
  }
  return { x: minX, y: minY, i: ind }
}
