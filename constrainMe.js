import Delaunator from 'delaunator'
import { heapSort } from './helpers'

// const test = [10, 6, 3, 4, 7, 1, 2, 5]

class Boundary {
  constructor (arr) {
    this.coords = arr.slice()
    this.index = [...this.coords.keys()].filter((i) => i % 2 === 0)
    this.clean()
    this.center = this.calcCenter()
    this.minY = minimumPointY(this.coords, this.index)
    this.sortHeapAndClean(this.coords, this.index, 'polar', [this.center.x, this.center.y])
    this.minY = minimumPointY(this.coords, this.index)
    this.minX = minimumPointX(this.coords, this.index)
    this.ray = null
    this.hull = []

    // center of test from html is not inside boundary
    // this point is though
    this.testFunctions()
  }

  testFunctions () {
    this.pointInOrOut([this.center.x, this.center.y])
    this.pointInOrOut([this.minX.x + 1000, this.minX.y])
    console.log(this.pointInOrOut([180, 100]))
   // this.concave(10)
  }

  concave (k) {
    // k nearest neighbor babbbbyyyy
    // https://towardsdatascience.com/the-concave-hull-c649795c0f0f
    // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
    // double check arr is sorted and clean
    this.sortHeapAndClean(this.coords, 2, 'polar', [this.center.x, this.center.y])

    let index = this.index.slice()
    if (index.length < 3) {
      return null
    } else if (index.length === 3) {
      return index
    }
    const kk = Math.min(Math.max(k, 3), index.length - 1)
    const firstPoint = this.minY.i
    const currentPoint = firstPoint
    const previousAngle = 0
    this.hull = [firstPoint]

    // TODO why is step init to 2?
    const step = 2
    index = index.splice(index.indexOf(firstPoint), 1)

    while ((currentPoint !== firstPoint || step === 2) && (index.length > 0)) {
      if (step === 5) {
        index = index.push(firstPoint)
      }
      const kNearestPoints = this.nearestPoints(index, currentPoint, kk)
    }
  }

  nearestPoints (index, cP, kk) {
    
  }

  sortHeapAndClean (arr, ind, criteria, centerPoint) {
    // console.log(this.index, this.coords2D)
    this.index = sortHeap(arr.slice(), this.index.slice(), criteria, centerPoint)
    // console.log(this.index, this.coords2D)
    this.clean()
    return this.coords
  }

  clean () {
    // TODO there has to be a better way to do this
    // On^2  urrrgh
    const index = this.index.slice()

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
    this.index = newIndex
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

  pointInOrOut (point) {
    // assume ray going to + infinity on x plane here just making assumption that it extends 1000 units past whatever the minimum x value is in the boundary
    const p = {
      x0: point[0], y0: point[1], x1: this.minX.x - 1000, y1: point[1]
    }
    this.ray = p
    // console.log(this.ray)
    // lets use non-zero winding number rule
    let windingNum = 0

    for (let i = 0; i < this.index.length; i++) {
      const l = {
        x0: this.coords[this.index[i]],
        y0: this.coords[this.index[i] + 1],
        x1: this.coords[this.index[(i + 1) > this.index.length - 1 ? 0 : i + 1]],
        y1: this.coords[this.index[(i + 1) > this.index.length - 1 ? 0 : i + 1] + 1]
      }
      const intersect = this.intersect(p, l)
      if (isFinite(intersect.x)) {
        if (l.y1 - l.y0 > 0) {
          windingNum++
        } else if (l.y1 - l.y0 < 0) {
          windingNum--
        }
      }
    }
    // console.log(windingNum)
    return windingNum !== 0
  }

  intersect (p, l) {
    const den = ((l.y1 - l.y0) * (p.x1 - p.x0)) - ((l.x1 - l.x0) * (p.y1 - p.y0))
    if (den === 0) {
      return { x: Infinity, y: Infinity }
    }

    let a = p.y0 - l.y0
    let b = p.x0 - l.x0

    const num1 = ((l.x1 - l.x0) * a) - ((l.y1 - l.y0) * b)
    const num2 = ((p.x1 - p.x0) * a) - ((p.y1 - p.y0) * b)

    a = num1 / den
    b = num2 / den

    const rv = {
      x: p.x0 + (a * (p.x1 - p.x0)),
      y: p.y0 + (a * (p.y1 - p.y0))
    }

    const t = { a: false, b: false }

    // if (p.y1 === rv.y) {
    // console.log(a, b);
    // }
    //
    if (a >= 0 && a < 1) {
      t.a = true
    }
    if (b >= 0 && b < 1) {
      t.b = true
    }
    if (t.a && t.b) {
      return rv
    }
    return { x: Infinity, y: Infinity }
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
      // sortHeap(test, 1)
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

function sortHeap (arr, index, criteria, centerPoint) {
  // convert point arr to 2d -> easier for me to get my head around sorting

  const minPoint = minimumPointY(arr, index)
  // builtInSort([minX, minY], newArr);
  heapSort([minPoint.x, minPoint.y], index, arr, index.length, criteria, centerPoint)

  return index
}

function minimumPointY (newArr, index) {
  let ind = 0
  let minY = Infinity
  let minX = Infinity
  if (index) {
    for (let p = 0; p < index.length; p++) {
      if (newArr[p + 1] < minY) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = p
      } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = p
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

function minimumPointX (newArr, index) {
  let ind = 0
  let minY = Infinity
  let minX = Infinity
  if (index) {
    for (let p = 0; p < index.length; p++) {
      if (newArr[p] < minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = p
      } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
        minX = newArr[p]
        minY = newArr[p + 1]
        ind = p
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
