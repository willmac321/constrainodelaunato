var c

function nextHalfEdge (e) {
  return (e % 3 === 2) ? e - 2 : e + 1
}

function dotProduct (a, b) {
  const p = { x: a[0], y: a[1] }
  const o = { x: b[0], y: b[1] }
  return p.x * o.x + p.y * o.y
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

function builtInSort (point, arr) {
  return arr.sort((a, b) => {
    // return manhattenDist(point, a) - manhattenDist(point, b)
    return (dotProduct(point, a) - dotProduct(point, b))
  })
}

// heap sort 2d array by angle
export function heapSort (point, index, a, count, p, center) {
  let func
  if (p === 'dist') {
    func = manhattenDist
  } else if (p === 'polar') {
    func = dotPolar
    c = { x: center[0], y: center[1] }
  } else if (p === 'dot') {
    func = dotProduct
  } else if (!Array.isArray(a[0])) {
    point = point[0]
    func = (a, b) => a - b
  } else {
    func = dotProduct
  }

  heapify(point, a, index, count, func)

  let end = count - 1
  while (end > 0) {
    swap(index, end, 0)
    end--
    siftDown(point, a, index, 0, end, func)
  }
}

function heapify (point, a, index, count, func) {
  const par = (i) => Math.floor((i - 1) / 2)
  let start = par(count - 1)
  while (start >= 0) {
    siftDown(point, a, index, start, count - 1, func)
    start--
  }
}

function siftDown (point, a, index, start, end, func) {
  let root = start
  const left = (i) => 2 * i + 1

  while (left(root) <= end) {
    const child = left(root)
    let s = root
    const dot = (i) => {
      const t = func([a[index[i]], a[index[i] + 1]], point)
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
