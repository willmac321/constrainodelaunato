import { dotProduct, heapSort, manhattenDist, intersect, slope, swap } from './helpers'

var counter = 0

export default class Boundary {
  constructor (arr, k = 3) {
    this.k = k
    this.coords = arr.slice()
    this.index = [...this.coords.keys()].filter((i) => i % 2 === 0)
    this.index = this.clean(this.index)
    this.center = this.calcCenter()
    this.minY = minimumPointY(this.coords, this.index)
    this.minX = minimumPointX(this.coords, this.index)
    this.maxX = maximumPointX(this.coords, this.index)

    this.ray = null
    this.hull = this.findConcaveHull(k)
    this.index = this.hull
  }

  testFunctions () {
    this.pointInOrOut([this.center.x, this.center.y], this.index, this.minX.x - 10)
    this.pointInOrOut([this.minX.x + 1000, this.minX.y], this.index, this.minX.x - 10)
    this.index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minY.x, this.minY.y], [this.minX.x, this.minY.y])
  }

  findConcaveHull (k) {
    // alt index is sorted to minX value
    this.index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y])
    this.hull = this.concave(this.index.slice(), k)
    return this.hull
  }

  concave (index, k) {
    // k nearest neighbor babbbbyyyy
    // https://towardsdatascience.com/the-concave-hull-c649795c0f0f
    // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
    // double check arr is sorted and clean
    // also sort it so all points are in order from some min point  on the xy plane
    const stopVal = 192 // Infinity // 76 // Infinity // and beyond
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
      let its = true
      let i = -1
      while (its && i < cPoints.length - 1) {
        // This is so that when the first point is added to the end of the hull, it doesn't get used to check for intersections
        let lastPoint = 0
        if (cPoints[i] === firstPoint.coord) {
          lastPoint = 1
        }
        let j = 1
        its = false
        while (!its && j < hull.length - lastPoint) {
          const l = {
            x0: this.coords[hull[step - 1]],
            y0: this.coords[hull[step - 1] + 1],
            x1: this.coords[cPoints[i + 1]],
            y1: this.coords[cPoints[i + 1] + 1],
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
          // (p.x0 !== l.x0 && p.y0 !== l.y0) ||
          // if (l.x0 === 216 && l.y0 === 133) {
          //   console.log(l, p, ints, isFinite(ints.x), (p.x1 === l.x0 && p.y1 === l.y0))
          // }
          if (isFinite(ints.x) && !endpointsMatch) {
            console.log(l, p, ints, isFinite(ints.x), !endpointsMatch)
            its = true
          }
          j++
        }
        i++
      }
      if (its) {
        console.log('intersection found at k ', k, its)
        return this.concave(oldIndex, ++kk)
      }
      currentPoint = cPoints[i]
      hull.push(currentPoint)
      if (counter > stopVal) {
        return hull// .concat(cPoints)
      }
      index.splice(index.indexOf(currentPoint), 1)
      step++
    }
    let allInside = true
    for (const i of index) {
      allInside = this.pointInOrOut(
        [this.coords[i], this.coords[i + 1]],
        hull, this.maxX.x + 10)
      if (!allInside) {
        break
      }
    }
    if (!allInside) {
      console.log('Another time round')
      return this.concave(oldIndex, ++kk)
    }
    console.log('made it out')
    this.k = kk
    return hull
  }

  sortByAngle (kNearestPoints, currentPoint, lastPoint) {
    const lastPointIndex = lastPoint
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
    //    console.log(`current slope ${slope(currentPointArr, [this.coords[rv[0]], this.coords[rv[0] + 1]])} for ${currentPointArr} and ${[this.coords[rv[0]], this.coords[rv[0] + 1]]}`)

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
      const newDist = manhattenDist(currentPointArr, newPoint)
      // console.log(`point ${k} at slope ${slope(lastPoint, newPoint)} for ${lastPoint} and ${newPoint}`)
      // console.log(`new point ${dotProduct(lastPoint, newPoint)}`)
      if (lastSlope && lastDist && Math.abs(newSlope) === Math.abs(lastSlope) && newDist < lastDist) {
        // flipflop the two points in array order if the slopes are the same
        // sort by euclid instead of straight swap
        swap(rv, k, k - 1)
        lastDist = manhattenDist(currentPointArr, [this.coords[rv[k]], this.coords[rv[k] + 1]])
        // if ((this.ray.x1 === 220 && this.ray.y1 === 89) || (this.ray.x0 === 220 && this.ray.y0 === 89)) {
        //   console.log(this.ray, this.subset(rv))
        // }
      } else {
        lastDist = newDist
      }

      lastSlope = slope(lastPoint, newPoint)
    }

    return rv // sortHeap(this.coords, kNearestPoints, 'polar', lastPoint, currentPointArr )
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
  // console.log({ x: minX, y: minY, i: ind })
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
