var c

// Unused functions
/* eslint-disable */
function nextHalfEdge (e) {
  return (e % 3 === 2) ? e - 2 : e + 1
}

function builtInSort (point, arr) {
  return arr.sort((a, b) => {
    // return manhattenDist(point, a) - manhattenDist(point, b)
    return (dotProduct(point, a) - dotProduct(point, b))
  })
}

/* eslint-enable */

/***
 * intersect compares two lines, does not include endpoints!
 * @param p is a line segment of type {x0, y0, x1, y1}
 * @param l is a line segment of type {x0, y0, x1, y1}
 * @param checkEndpoints bool - True to include endpoints in intersection calc
 * @return {x, y} at intersect, if they dont intersect, they intersect at infinity
 ***/
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
  if (a >= 0 && a < 1) {
    t.a = true
  }
  if (b >= 0 && b < 1) {
    t.b = true
  }
  return t
}

function dotPolar (a, b) {
  // b is basis point
  // put those bad boys in order ccw around some centroid point that globally declared
  const o = { x: c.x - a[0], y: c.y - a[1] }
  const p = { x: c.x - b[0], y: c.y - b[1] }
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

function euclid (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return Math.sqrt(Math.pow(p.x - o.x, 2) + Math.pow(p.y - o.y, 2))
}

// heap sort 2d array by angle
export function heapSort (minpoint, index, a, count, p, center) {
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
//    console.log(i, a[i], a[i + 1], func([a[i], a[i + 1]], minpoint))
//  }
//  console.log(minpoint)
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

function swap (a, i, j) {
  const t = a[i]
  a[i] = a[j]
  a[j] = t
}
