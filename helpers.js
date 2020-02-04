var c

// Unused functions
/* eslint-disable */
function builtInSort (point, arr) {
  return arr.sort((a, b) => {
    // return manhattenDist(point, a) - manhattenDist(point, b)
    return (dotProduct(point, a) - dotProduct(point, b))
  })
}

function intersectAlt (p, l) {
  // http://bl.ocks.org/nitaku/fdbb70c3baa36e8feb4e
  const s1_x = p.x1 - p.x0
  const s1_y = p.y1 - p.y0
  const s2_x = l.x1 - l.x0
  const s2_y = l.y1 - l.y0

  const s = (-s1_y * (p.x0 - l.x0) + s1_x * (p.y0 - l.y0)) / (-s2_x * s1_y + s1_x * s2_y)
  const t = (s2_x * (p.y0 - l.y0) - s2_y * (p.x0 - l.x0)) / (-s2_x * s1_y + s1_x * s2_y)

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    return { x: p.x0 + (t * s1_x), y: p.y0 + (t * s1_y) }
  }
  return { x: Infinity, y: Infinity }
}

function polarAngle (a, b) {
  const o = { x: a[0] - c.x, y: a[1] - c.y }
  const p = { x: b[0] - c.x, y: b[1] - c.y }
  let theta = (Math.atan2(p.y, p.x) - (Math.atan2(o.y, o.x)) * 180 / Math.PI) % 360
  theta = isNaN(theta) ? 0 : theta
  theta = theta > 0 ? theta : 360 + theta
  return theta
}


/* eslint-enable */

function nextHalfEdge (e) {
  return (e % 3 === 2) ? e - 2 : e + 1
}

/**
 * getEdges
 *
 * @param {Object} delaunay Delaunator object
 * @returns {Array} array of indices for edge points, so the indices for a coord array that are in order of triangulation
 */
export function getEdges (delaunay) {
  const rv = []
  for (let e = 0; e < delaunay.triangles.length; e++) {
    if (e > delaunay.halfedges[e]) {
      rv.push(2 * delaunay.triangles[e], 2 * delaunay.triangles[nextHalfEdge(e)])
    }
  }
  return rv
}

/**
 * intersect
 * compares two lines, does not include endpoints!
 * @param {Object} p object with 4 int values of type {x0, y0, x1, y1} where 0 denotes start point and 1 is endpoint
 * @param {Object} l object with 4 int values of type {x0, y0, x1, y1} where 0 denotes start point and 1 is endpoint
 * @param {bool} checkEndpoints=false whether or not to check endpoints in intersection calc, does not check endpoints by default, true to check them
 * @returns {Object} Intersection point x and y coords, returns x: Inf and y: Inf if points do not intersect
 */
export function intersect (p, l, checkEndpoints = false) {
  // compare two line segments to see if they intersect
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

  // if (p.y1 === rv.y) {
  // console.log(a, b);
  // }
  //
  let t = compareIntersect(a, b)
  if (checkEndpoints) {
    t = compareIntersectEndpoints(a, b)
  }

  if (t.a && t.b) {
    return rv
  }
  return { x: Infinity, y: Infinity }
}

export function dotProduct (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return p.x * o.x + p.y * o.y
}

export function slope (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return (p.y - o.y) / (p.x - o.x)
}

function compareIntersect (a, b) {
  const t = { a: false, b: false }
  if (a > 0 && a < 1) {
    t.a = true
  }
  if (b > 0 && b < 1) {
    t.b = true
  }
  return t
}

function compareIntersectEndpoints (a, b) {
  const t = { a: false, b: false }
  if (a >= 0 && a <= 1) {
    t.a = true
  }
  if (b >= 0 && b <= 1) {
    t.b = true
  }
  return t
}

function dotPolar (a, b) {
  // b is basis point
  // put those bad boys in order ccw around some centroid point that globally declared
  const o = { x: a[0] - c.x, y: a[1] - c.y }
  const p = { x: b[0] - c.x, y: b[1] - c.y }
  const magO = Math.sqrt(Math.pow(o.x, 2) + Math.pow(o.y, 2))
  const magP = Math.sqrt(Math.pow(p.x, 2) + Math.pow(p.y, 2))
  let theta = Math.acos((p.x * o.x + p.y * o.y) / (magO * magP)) * 180 / Math.PI
  const det = (o.x * p.y) - (p.x * o.y)
  theta = isNaN(theta) ? 0 : theta
  theta = det > 0 ? 360 - theta : theta
  return theta
}

function manhattenDist (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return Math.abs(p.x - o.x) + Math.abs(p.y - o.y)
}

/**
 * euclid
 *
 * @param {Array} a x and y coord with a[0] as x
 * @param {Array} b x and y coord with b[0] as x
 * @returns {Double} Float/Double value of euclidian distance between two points
 */
export function euclid (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return Math.sqrt(Math.pow(p.x - o.x, 2) + Math.pow(p.y - o.y, 2))
}

/**
 * distLineAndPoint
 * calculates the distance between a point and a line derived from a supplied line segment
 * this distance may intersect the line outside of the actual line segment
 * https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
 *
 * @param {Object} l Line Segment format {x0, y0, x1, y1}
 * @param {Object} p point format {x, y}
 * @returns {Float} distance between point and line
 */
export function distLineAndPoint (l, p) {
  const num = Math.abs((l.y1 - l.y0) * p.x - (l.x1 - l.x0) * p.y + l.x1 * l.y0 - l.y1 * l.x0)
  const dem = Math.sqrt(Math.pow(l.y1 - l.y0, 2) + Math.pow(l.x1 - l.x0, 2))
  return num / dem
}

// heap sort 2d array by angle
function heapSort (minpoint, index, a, count, p, center) {
  let func

  if (p === 'dist') {
    func = manhattenDist
  } else if (p === 'distrel') {
    func = manhattenDist
  } else if (p === 'euclid') {
    func = euclid
  } else if (p === 'polar') {
    func = dotPolar
    c = { x: center[0], y: center[1] }
  } else if (p === 'dot') {
    func = dotProduct
  } else if (!Array.isArray(a[0])) {
    minpoint = minpoint[0]
    func = (a, b) => a - b
  } else {
    func = dotProduct
  }

  heapify(minpoint, a, index, count, func, p)

  let end = count - 1
  while (end > 0) {
    swap(index, end, 0)
    end--
    siftDown(minpoint, a, index, 0, end, func, p)
  }
//  for (const i of index) {
//    console.log(i, a[i], a[i + 1], func([a[i], a[i + 1]], minpoint), minpoint)
//  }
}

function heapify (point, a, index, count, func, p) {
  const par = (i) => Math.floor((i - 1) / 2)
  let start = par(count - 1)
  while (start >= 0) {
    siftDown(point, a, index, start, count - 1, func, p)
    start--
  }
}

function siftDown (point, a, index, start, end, func, p) {
  let root = start
  const left = (i) => 2 * i + 1

  while (left(root) <= end) {
    const child = left(root)
    let s = root
    const dot = (i) => {
      const t = func([a[index[i]], a[index[i] + 1]], point)
      if (p === 'distrel') {
        point = [a[index[i]], a[index[i] + 1]]
      }
      return t
    }

    if (dot(s) < dot(child)) {
      s = child
    }
    if (child + 1 <= end && dot(s) < dot(child + 1)) {
      s = child + 1
    }
    if (s === root) {
      return
    } else {
      swap(index, root, s)
      root = s
    }
  }
}

export function swap (a, i, j) {
  const t = a[i]
  a[i] = a[j]
  a[j] = t
}

export function maximumPointX (newArr, index) {
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

export function minimumPointY (newArr, index) {
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

export function sortHeap (arr, index, criteria, minPoint, centerPoint) {
  // convert point arr to 2d -> easier for me to get my head around sorting

  // minPoint = { x: minPoint.x, y: minPoint.y }
  // builtInSort([minX, minY], newArr);

  heapSort(minPoint, index, arr, index.length, criteria, centerPoint)

  return index
}

export function minimumPointX (newArr, index) {
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
