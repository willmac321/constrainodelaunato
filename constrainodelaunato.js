(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global = global || self, global.ConstrainoDelaunato = factory());
}(this, (function () { 'use strict';

    const EPSILON = Math.pow(2, -52);
    const EDGE_STACK = new Uint32Array(512);

    class Delaunator {

        static from(points, getX = defaultGetX, getY = defaultGetY) {
            const n = points.length;
            const coords = new Float64Array(n * 2);

            for (let i = 0; i < n; i++) {
                const p = points[i];
                coords[2 * i] = getX(p);
                coords[2 * i + 1] = getY(p);
            }

            return new Delaunator(coords);
        }

        constructor(coords) {
            const n = coords.length >> 1;
            if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.');

            this.coords = coords;

            // arrays that will store the triangulation graph
            const maxTriangles = Math.max(2 * n - 5, 0);
            this._triangles = new Uint32Array(maxTriangles * 3);
            this._halfedges = new Int32Array(maxTriangles * 3);

            // temporary arrays for tracking the edges of the advancing convex hull
            this._hashSize = Math.ceil(Math.sqrt(n));
            this._hullPrev = new Uint32Array(n); // edge to prev edge
            this._hullNext = new Uint32Array(n); // edge to next edge
            this._hullTri = new Uint32Array(n); // edge to adjacent triangle
            this._hullHash = new Int32Array(this._hashSize).fill(-1); // angular edge hash

            // temporary arrays for sorting points
            this._ids = new Uint32Array(n);
            this._dists = new Float64Array(n);

            this.update();
        }

        update() {
            const {coords, _hullPrev: hullPrev, _hullNext: hullNext, _hullTri: hullTri, _hullHash: hullHash} =  this;
            const n = coords.length >> 1;

            // populate an array of point indices; calculate input data bbox
            let minX = Infinity;
            let minY = Infinity;
            let maxX = -Infinity;
            let maxY = -Infinity;

            for (let i = 0; i < n; i++) {
                const x = coords[2 * i];
                const y = coords[2 * i + 1];
                if (x < minX) minX = x;
                if (y < minY) minY = y;
                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
                this._ids[i] = i;
            }
            const cx = (minX + maxX) / 2;
            const cy = (minY + maxY) / 2;

            let minDist = Infinity;
            let i0, i1, i2;

            // pick a seed point close to the center
            for (let i = 0; i < n; i++) {
                const d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist) {
                    i0 = i;
                    minDist = d;
                }
            }
            const i0x = coords[2 * i0];
            const i0y = coords[2 * i0 + 1];

            minDist = Infinity;

            // find the point closest to the seed
            for (let i = 0; i < n; i++) {
                if (i === i0) continue;
                const d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
                if (d < minDist && d > 0) {
                    i1 = i;
                    minDist = d;
                }
            }
            let i1x = coords[2 * i1];
            let i1y = coords[2 * i1 + 1];

            let minRadius = Infinity;

            // find the third point which forms the smallest circumcircle with the first two
            for (let i = 0; i < n; i++) {
                if (i === i0 || i === i1) continue;
                const r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
                if (r < minRadius) {
                    i2 = i;
                    minRadius = r;
                }
            }
            let i2x = coords[2 * i2];
            let i2y = coords[2 * i2 + 1];

            if (minRadius === Infinity) {
                // order collinear points by dx (or dy if all x are identical)
                // and return the list as a hull
                for (let i = 0; i < n; i++) {
                    this._dists[i] = (coords[2 * i] - coords[0]) || (coords[2 * i + 1] - coords[1]);
                }
                quicksort(this._ids, this._dists, 0, n - 1);
                const hull = new Uint32Array(n);
                let j = 0;
                for (let i = 0, d0 = -Infinity; i < n; i++) {
                    const id = this._ids[i];
                    if (this._dists[id] > d0) {
                        hull[j++] = id;
                        d0 = this._dists[id];
                    }
                }
                this.hull = hull.subarray(0, j);
                this.triangles = new Uint32Array(0);
                this.halfedges = new Uint32Array(0);
                return;
            }

            // swap the order of the seed points for counter-clockwise orientation
            if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
                const i = i1;
                const x = i1x;
                const y = i1y;
                i1 = i2;
                i1x = i2x;
                i1y = i2y;
                i2 = i;
                i2x = x;
                i2y = y;
            }

            const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
            this._cx = center.x;
            this._cy = center.y;

            for (let i = 0; i < n; i++) {
                this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
            }

            // sort the points by distance from the seed triangle circumcenter
            quicksort(this._ids, this._dists, 0, n - 1);

            // set up the seed triangle as the starting hull
            this._hullStart = i0;
            let hullSize = 3;

            hullNext[i0] = hullPrev[i2] = i1;
            hullNext[i1] = hullPrev[i0] = i2;
            hullNext[i2] = hullPrev[i1] = i0;

            hullTri[i0] = 0;
            hullTri[i1] = 1;
            hullTri[i2] = 2;

            hullHash.fill(-1);
            hullHash[this._hashKey(i0x, i0y)] = i0;
            hullHash[this._hashKey(i1x, i1y)] = i1;
            hullHash[this._hashKey(i2x, i2y)] = i2;

            this.trianglesLen = 0;
            this._addTriangle(i0, i1, i2, -1, -1, -1);

            for (let k = 0, xp, yp; k < this._ids.length; k++) {
                const i = this._ids[k];
                const x = coords[2 * i];
                const y = coords[2 * i + 1];

                // skip near-duplicate points
                if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
                xp = x;
                yp = y;

                // skip seed triangle points
                if (i === i0 || i === i1 || i === i2) continue;

                // find a visible edge on the convex hull using edge hash
                let start = 0;
                for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
                    start = hullHash[(key + j) % this._hashSize];
                    if (start !== -1 && start !== hullNext[start]) break;
                }

                start = hullPrev[start];
                let e = start, q;
                while (q = hullNext[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
                    e = q;
                    if (e === start) {
                        e = -1;
                        break;
                    }
                }
                if (e === -1) continue; // likely a near-duplicate point; skip it

                // add the first triangle from the point
                let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

                // recursively flip triangles from the point until they satisfy the Delaunay condition
                hullTri[i] = this._legalize(t + 2);
                hullTri[e] = t; // keep track of boundary triangles on the hull
                hullSize++;

                // walk forward through the hull, adding more triangles and flipping recursively
                let n = hullNext[e];
                while (q = hullNext[n], orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
                    t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                    hullTri[i] = this._legalize(t + 2);
                    hullNext[n] = n; // mark as removed
                    hullSize--;
                    n = q;
                }

                // walk backward from the other side, adding more triangles and flipping
                if (e === start) {
                    while (q = hullPrev[e], orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
                        t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                        this._legalize(t + 2);
                        hullTri[q] = t;
                        hullNext[e] = e; // mark as removed
                        hullSize--;
                        e = q;
                    }
                }

                // update the hull indices
                this._hullStart = hullPrev[i] = e;
                hullNext[e] = hullPrev[n] = i;
                hullNext[i] = n;

                // save the two new edges in the hash table
                hullHash[this._hashKey(x, y)] = i;
                hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
            }

            this.hull = new Uint32Array(hullSize);
            for (let i = 0, e = this._hullStart; i < hullSize; i++) {
                this.hull[i] = e;
                e = hullNext[e];
            }

            // trim typed triangle mesh arrays
            this.triangles = this._triangles.subarray(0, this.trianglesLen);
            this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
        }

        _hashKey(x, y) {
            return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
        }

        _legalize(a) {
            const {_triangles: triangles, _halfedges: halfedges, coords} = this;

            let i = 0;
            let ar = 0;

            // recursion eliminated with a fixed-size stack
            while (true) {
                const b = halfedges[a];

                /* if the pair of triangles doesn't satisfy the Delaunay condition
                 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
                 * then do the same check/flip recursively for the new pair of triangles
                 *
                 *           pl                    pl
                 *          /||\                  /  \
                 *       al/ || \bl            al/    \a
                 *        /  ||  \              /      \
                 *       /  a||b  \    flip    /___ar___\
                 *     p0\   ||   /p1   =>   p0\---bl---/p1
                 *        \  ||  /              \      /
                 *       ar\ || /br             b\    /br
                 *          \||/                  \  /
                 *           pr                    pr
                 */
                const a0 = a - a % 3;
                ar = a0 + (a + 2) % 3;

                if (b === -1) { // convex hull edge
                    if (i === 0) break;
                    a = EDGE_STACK[--i];
                    continue;
                }

                const b0 = b - b % 3;
                const al = a0 + (a + 1) % 3;
                const bl = b0 + (b + 2) % 3;

                const p0 = triangles[ar];
                const pr = triangles[a];
                const pl = triangles[al];
                const p1 = triangles[bl];

                const illegal = inCircle(
                    coords[2 * p0], coords[2 * p0 + 1],
                    coords[2 * pr], coords[2 * pr + 1],
                    coords[2 * pl], coords[2 * pl + 1],
                    coords[2 * p1], coords[2 * p1 + 1]);

                if (illegal) {
                    triangles[a] = p1;
                    triangles[b] = p0;

                    const hbl = halfedges[bl];

                    // edge swapped on the other side of the hull (rare); fix the halfedge reference
                    if (hbl === -1) {
                        let e = this._hullStart;
                        do {
                            if (this._hullTri[e] === bl) {
                                this._hullTri[e] = a;
                                break;
                            }
                            e = this._hullPrev[e];
                        } while (e !== this._hullStart);
                    }
                    this._link(a, hbl);
                    this._link(b, halfedges[ar]);
                    this._link(ar, bl);

                    const br = b0 + (b + 1) % 3;

                    // don't worry about hitting the cap: it can only happen on extremely degenerate input
                    if (i < EDGE_STACK.length) {
                        EDGE_STACK[i++] = br;
                    }
                } else {
                    if (i === 0) break;
                    a = EDGE_STACK[--i];
                }
            }

            return ar;
        }

        _link(a, b) {
            this._halfedges[a] = b;
            if (b !== -1) this._halfedges[b] = a;
        }

        // add a new triangle given vertex indices and adjacent half-edge ids
        _addTriangle(i0, i1, i2, a, b, c) {
            const t = this.trianglesLen;

            this._triangles[t] = i0;
            this._triangles[t + 1] = i1;
            this._triangles[t + 2] = i2;

            this._link(t, a);
            this._link(t + 1, b);
            this._link(t + 2, c);

            this.trianglesLen += 3;

            return t;
        }
    }

    // monotonically increases with real angle, but doesn't need expensive trigonometry
    function pseudoAngle(dx, dy) {
        const p = dx / (Math.abs(dx) + Math.abs(dy));
        return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
    }

    function dist(ax, ay, bx, by) {
        const dx = ax - bx;
        const dy = ay - by;
        return dx * dx + dy * dy;
    }

    // return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
    function orientIfSure(px, py, rx, ry, qx, qy) {
        const l = (ry - py) * (qx - px);
        const r = (rx - px) * (qy - py);
        return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
    }

    // a more robust orientation test that's stable in a given triangle (to fix robustness issues)
    function orient(rx, ry, qx, qy, px, py) {
        const sign = orientIfSure(px, py, rx, ry, qx, qy) ||
        orientIfSure(rx, ry, qx, qy, px, py) ||
        orientIfSure(qx, qy, px, py, rx, ry);
        return sign < 0;
    }

    function inCircle(ax, ay, bx, by, cx, cy, px, py) {
        const dx = ax - px;
        const dy = ay - py;
        const ex = bx - px;
        const ey = by - py;
        const fx = cx - px;
        const fy = cy - py;

        const ap = dx * dx + dy * dy;
        const bp = ex * ex + ey * ey;
        const cp = fx * fx + fy * fy;

        return dx * (ey * cp - bp * fy) -
               dy * (ex * cp - bp * fx) +
               ap * (ex * fy - ey * fx) < 0;
    }

    function circumradius(ax, ay, bx, by, cx, cy) {
        const dx = bx - ax;
        const dy = by - ay;
        const ex = cx - ax;
        const ey = cy - ay;

        const bl = dx * dx + dy * dy;
        const cl = ex * ex + ey * ey;
        const d = 0.5 / (dx * ey - dy * ex);

        const x = (ey * bl - dy * cl) * d;
        const y = (dx * cl - ex * bl) * d;

        return x * x + y * y;
    }

    function circumcenter(ax, ay, bx, by, cx, cy) {
        const dx = bx - ax;
        const dy = by - ay;
        const ex = cx - ax;
        const ey = cy - ay;

        const bl = dx * dx + dy * dy;
        const cl = ex * ex + ey * ey;
        const d = 0.5 / (dx * ey - dy * ex);

        const x = ax + (ey * bl - dy * cl) * d;
        const y = ay + (dx * cl - ex * bl) * d;

        return {x, y};
    }

    function quicksort(ids, dists, left, right) {
        if (right - left <= 20) {
            for (let i = left + 1; i <= right; i++) {
                const temp = ids[i];
                const tempDist = dists[temp];
                let j = i - 1;
                while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];
                ids[j + 1] = temp;
            }
        } else {
            const median = (left + right) >> 1;
            let i = left + 1;
            let j = right;
            swap(ids, median, i);
            if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
            if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
            if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

            const temp = ids[i];
            const tempDist = dists[temp];
            while (true) {
                do i++; while (dists[ids[i]] < tempDist);
                do j--; while (dists[ids[j]] > tempDist);
                if (j < i) break;
                swap(ids, i, j);
            }
            ids[left + 1] = ids[j];
            ids[j] = temp;

            if (right - i + 1 >= j - left) {
                quicksort(ids, dists, i, right);
                quicksort(ids, dists, left, j - 1);
            } else {
                quicksort(ids, dists, left, j - 1);
                quicksort(ids, dists, i, right);
            }
        }
    }

    function swap(arr, i, j) {
        const tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    function defaultGetX(p) {
        return p[0];
    }
    function defaultGetY(p) {
        return p[1];
    }

    var c;


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
    function getEdges (delaunay) {
      const rv = [];
      for (let e = 0; e < delaunay.triangles.length; e++) {
        if (e > delaunay.halfedges[e]) {
          rv.push(2 * delaunay.triangles[e], 2 * delaunay.triangles[nextHalfEdge(e)]);
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
    function intersect (p, l, checkEndpoints = false) {
      // compare two line segments to see if they intersect
      const den = ((l.y1 - l.y0) * (p.x1 - p.x0)) - ((l.x1 - l.x0) * (p.y1 - p.y0));
      if (den === 0) {
        return { x: Infinity, y: Infinity }
      }

      let a = p.y0 - l.y0;
      let b = p.x0 - l.x0;

      const num1 = ((l.x1 - l.x0) * a) - ((l.y1 - l.y0) * b);
      const num2 = ((p.x1 - p.x0) * a) - ((p.y1 - p.y0) * b);

      a = num1 / den;
      b = num2 / den;

      const rv = {
        x: p.x0 + (a * (p.x1 - p.x0)),
        y: p.y0 + (a * (p.y1 - p.y0))
      };

      // if (p.y1 === rv.y) {
      // console.log(a, b);
      // }
      //
      let t = compareIntersect(a, b);
      if (checkEndpoints) {
        t = compareIntersectEndpoints(a, b);
      }

      if (t.a && t.b) {
        return rv
      }
      return { x: Infinity, y: Infinity }
    }

    function dotProduct (a, b) {
      const p = { x: a[0], y: a[1] };
      const o = { x: b[0], y: b[1] };
      return p.x * o.x + p.y * o.y
    }

    function slope (a, b) {
      const p = { x: a[0], y: a[1] };
      const o = { x: b[0], y: b[1] };
      const rv = (p.y - o.y) / (p.x - o.x);

      return rv
    }

    function compareIntersect (a, b) {
      const t = { a: false, b: false };
      if (a > 0 && a < 1) {
        t.a = true;
      }
      if (b > 0 && b < 1) {
        t.b = true;
      }
      return t
    }

    function compareIntersectEndpoints (a, b) {
      const t = { a: false, b: false };
      if (a >= 0 && a <= 1) {
        t.a = true;
      }
      if (b >= 0 && b <= 1) {
        t.b = true;
      }
      return t
    }

    function dotPolar (a, b) {
      // b is basis point
      // put those bad boys in order ccw around some centroid point that globally declared
      const o = { x: a[0] - c.x, y: a[1] - c.y };
      const p = { x: b[0] - c.x, y: b[1] - c.y };
      const magO = Math.sqrt(Math.pow(o.x, 2) + Math.pow(o.y, 2));
      const magP = Math.sqrt(Math.pow(p.x, 2) + Math.pow(p.y, 2));
      let theta = Math.acos((p.x * o.x + p.y * o.y) / (magO * magP)) * 180 / Math.PI;
      const det = (o.x * p.y) - (p.x * o.y);
      theta = isNaN(theta) ? 0 : theta;
      theta = det > 0 ? 360 - theta : theta;
      return theta
    }

    function manhattenDist (a, b) {
      const p = { x: a[0], y: a[1] };
      const o = { x: b[0], y: b[1] };
      return Math.abs(p.x - o.x) + Math.abs(p.y - o.y)
    }

    /**
     * euclid
     *
     * @param {Array} a x and y coord with a[0] as x
     * @param {Array} b x and y coord with b[0] as x
     * @returns {Double} Float/Double value of euclidian distance between two points
     */
    function euclid (a, b) {
      const p = { x: a[0], y: a[1] };
      const o = { x: b[0], y: b[1] };
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
    function distLineAndPoint (l, p) {
      const num = Math.abs((l.y1 - l.y0) * p.x - (l.x1 - l.x0) * p.y + l.x1 * l.y0 - l.y1 * l.x0);
      const dem = Math.sqrt(Math.pow(l.y1 - l.y0, 2) + Math.pow(l.x1 - l.x0, 2));
      return num / dem
    }

    // heap sort 2d array by angle
    function heapSort (minpoint, index, a, count, p, center) {
      let func;

      if (p === 'dist') {
        func = manhattenDist;
      } else if (p === 'distrel') {
        func = manhattenDist;
      } else if (p === 'euclid') {
        func = euclid;
      } else if (p === 'polar') {
        func = dotPolar;
        c = { x: center[0], y: center[1] };
      } else if (p === 'dot') {
        func = dotProduct;
      } else if (!Array.isArray(a[0])) {
        minpoint = minpoint[0];
        func = (a, b) => a - b;
      } else {
        func = dotProduct;
      }

      heapify(minpoint, a, index, count, func, p);

      let end = count - 1;
      while (end > 0) {
        swap$1(index, end, 0);
        end--;
        siftDown(minpoint, a, index, 0, end, func, p);
      }
    //  for (const i of index) {
    //    console.log(i, a[i], a[i + 1], func([a[i], a[i + 1]], minpoint), minpoint)
    //  }
    }

    function heapify (point, a, index, count, func, p) {
      const par = (i) => Math.floor((i - 1) / 2);
      let start = par(count - 1);
      while (start >= 0) {
        siftDown(point, a, index, start, count - 1, func, p);
        start--;
      }
    }

    function siftDown (point, a, index, start, end, func, p) {
      let root = start;
      const left = (i) => 2 * i + 1;

      while (left(root) <= end) {
        const child = left(root);
        let s = root;
        const dot = (i) => {
          const t = func([a[index[i]], a[index[i] + 1]], point);
          if (p === 'distrel') {
            point = [a[index[i]], a[index[i] + 1]];
          }
          return t
        };

        if (dot(s) < dot(child)) {
          s = child;
        }
        if (child + 1 <= end && dot(s) < dot(child + 1)) {
          s = child + 1;
        }
        if (s === root) {
          return
        } else {
          swap$1(index, root, s);
          root = s;
        }
      }
    }

    function swap$1 (a, i, j) {
      const t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
    function maximumPointY (newArr, index) {
      let ind = 0;
      let minY = -Infinity;
      let minX = -Infinity;
      if (index) {
        for (const [k, p] of index.entries()) {
          if (newArr[p + 1] > minY) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          } else if (newArr[p + 1] >= minY && newArr[p] >= minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          }
        }
      } else {
        for (let p = 0; p < newArr.length; p++) {
          if (newArr[p] > minX) {
            minX = newArr[p];
            ind = p;
          }
        }
      }
      // console.log({ x: minX, y: minY, i: ind })
      return { x: minX, y: minY, i: ind }
    }

    function maximumPointX (newArr, index) {
      let ind = 0;
      let minY = -Infinity;
      let minX = -Infinity;
      if (index) {
        for (const [k, p] of index.entries()) {
          if (newArr[p] > minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          } else if (newArr[p + 1] >= minY && newArr[p] >= minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          }
        }
      } else {
        for (let p = 0; p < newArr.length; p++) {
          if (newArr[p] > minX) {
            minX = newArr[p];
            ind = p;
          }
        }
      }
      return { x: minX, y: minY, i: ind }
    }

    function minimumPointY (newArr, index) {
      let ind = 0;
      let minY = Infinity;
      let minX = Infinity;
      if (index) {
        for (const [k, p] of index.entries()) {
          if (newArr[p + 1] < minY) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          }
        }
      } else {
        for (let p = 0; p < newArr.length; p++) {
          if (newArr[p] < minX) {
            minX = newArr[p];
            ind = p;
          }
        }
      }
      // console.log({ x: minX, y: minY, i: ind })
      return { x: minX, y: minY, i: ind }
    }

    function sortHeap (arr, index, criteria, minPoint, centerPoint) {
      // convert point arr to 2d -> easier for me to get my head around sorting

      // minPoint = { x: minPoint.x, y: minPoint.y }
      // builtInSort([minX, minY], newArr);

      heapSort(minPoint, index, arr, index.length, criteria, centerPoint);

      return index
    }

    function minimumPointX (newArr, index) {
      let ind = 0;
      let minY = Infinity;
      let minX = Infinity;
      if (index) {
        for (const [k, p] of index.entries()) {
          if (newArr[p] < minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          } else if (newArr[p + 1] <= minY && newArr[p] <= minX) {
            minX = newArr[p];
            minY = newArr[p + 1];
            ind = k;
          }
        }
      } else {
        for (let p = 0; p < newArr.length; p++) {
          if (newArr[p] < minX) {
            minX = newArr[p];
            ind = p;
          }
        }
      }
      return { x: minX, y: minY, i: ind }
    }

    var counter = 0;

    class Boundary {
      constructor (arr, k = 3, dist) {
        this.k = k;
        this.dist = dist;
        this.coords = arr.slice();
        this.index = [...this.coords.keys()].filter((i) => i % 2 === 0);
        this.index = this.clean(this.index);
        this.center = this.calcCenter();
        this.minY = minimumPointY(this.coords, this.index);
        this.maxY = maximumPointY(this.coords, this.index);
        this.minX = minimumPointX(this.coords, this.index);
        this.maxX = maximumPointX(this.coords, this.index);
        this.maxD = Math.sqrt(Math.pow(this.maxX.x - this.minX.x, 2) + Math.pow(this.maxY.y - this.minY.y, 2));
        this.offsetAngle = 1;

        this.cPoints = [];

        this.ray = null;
        this.hull = this.findConcaveHull(k);
      }

      /**
       * findConcaveHull
       *
       * @param {Integer} k Starting point cloud count for points being used for selection
       * @returns {Array} Array of indices of X values for the sorted concave hull
       */
      findConcaveHull (k) {
        // alt index is sorted to minX value
        const index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y]);
        const indexCopy = [...index];
        let hull = this.concave(index, k);
        if (hull === null) {
          // try enlarging distance to find match
          this.dist = Infinity;
          hull = this.concave(indexCopy, k);
        }

        return hull
      }

      concave (index, k) {
        // k nearest neighbor babbbbyyyy
        // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
        // double check arr is sorted and clean
        // also sort it so all points are in order from some min point  on the xy plane
        const stopVal = Infinity; // 45 // Infinity // 200 // 86 // Infinity // 76 // Infinity // and beyond
        const oldIndex = index.slice();
        if (index.length < 3) {
          console.error(`Remaining points can not form a polygon k:${k}`);
          return null
        } else if (k > index.length - 1) {
          console.error(`K exceeds total length of points available k:${k}`);
          return null
        } else if (index.length === 3) {
          return index
        }

        let kk = Math.min(Math.max(k, 3), index.length - 1);
        // i is a pointer to the relative index not a loc in this.coords
        // so, index of that index gives a this.coords pointer
        const firstPointIndex = minimumPointY(this.coords, index).i;
        const firstPoint = { i: firstPointIndex, coord: index[firstPointIndex] };
        let currentPoint = firstPoint.coord;
        const hull = [firstPoint.coord];
        // why is step init to 2?
        // Because the paper was written in Matlab....
        let step = 1;
        // each index value can only be used once so this is ok
        index.splice(firstPoint.i, 1);
        while ((currentPoint !== firstPoint.coord || step === 1) && (index.length > 0)) {
          counter++;
          if (step === 4 && index.indexOf(firstPoint.coord) < 0) {
            index.push(firstPoint.coord);
          } else if (step <= oldIndex.length - 1 && step < 4 && (index.length === 1) && index.indexOf(firstPoint.coord) < 0) {
            // if the poly is 3 or 4 points add the index to the end so we dont short circuit
            index.push(firstPoint.coord);
          }
          // find nearest neighbors
          const kNearestPoints = this.nearestPoints(index, currentPoint, kk);
          // descending order 'right-hand' turn x and y min are top left on js canvas in webpage
          const cPoints = this.sortByAngle(kNearestPoints, currentPoint, hull[hull.length - 2]);
          // if (cPoints.indexOf(firstPoint.coord) > -1) {
          //   console.log(cPoints)
          // }
          let its = true;
          let i = -1;
          while (its && i < cPoints.length - 1) {
            // This is so that when the first point is added to the end of the hull, it doesn't get used to check for intersections
            let lastPoint = 0;
            if (cPoints[i] === firstPoint.coord) {
              // console.log('back to first', firstPoint)
              lastPoint = 1;
            }
            let j = 1;
            its = false;
            while (!its && j < hull.length - lastPoint) {
              const l = {
                x0: this.coords[hull[step - 1]],
                y0: this.coords[hull[step - 1] + 1],
                x1: this.coords[cPoints[i + 1]],
                y1: this.coords[cPoints[i + 1] + 1]
              };
              const p = {
                x0: this.coords[hull[step - j]],
                y0: this.coords[hull[step - j] + 1],
                x1: this.coords[hull[step - 1 - j]],
                y1: this.coords[hull[step - 1 - j] + 1]
              };
              // the endpoint of one line segment is always intersecting the endpoint of a connected line segment, how to ignore this intersection?
              const ints = intersect(p, l, true);
              const endpointsMatch = (p.x0 === l.x0 && p.y0 === l.y0);
              const isClose = (cPoints[i + 1] === firstPoint.coord) && (p.x1 === l.x1 && p.y1 === l.y1);
              // (p.x0 !== l.x0 && p.y0 !== l.y0) ||
              // if (l.x0 === 221 && l.y0 === 90) {
              //   console.log(l, p, ints, isFinite(ints.x), (p.x1 === l.x0 && p.y1 === l.y0))
              // }
              if (isFinite(ints.x) && !endpointsMatch && !isClose) {
                its = true;
              }
              j++;
            }
            i++;
          }
          this.cPoints = cPoints.slice();
          // this.cPoints.splice(i, 1)

          if (its) {
            return this.concave(oldIndex, ++kk)
          }

          currentPoint = cPoints[i];
          hull.push(currentPoint);

          if (counter > stopVal) {
            return hull // .concat(cPoints)
          }

          index.splice(index.indexOf(currentPoint), 1);

          step++;
        }

        let allInside = true;
        for (const i of index) {
          allInside = this.pointInOrOut(
            [this.coords[i], this.coords[i + 1]],
            hull, 0);
          if (!allInside) {
            break
          }
        }
        if (!allInside) {
          return this.concave(oldIndex, ++kk)
        }
        this.k = kk;
        return hull
      }

      sortByAngle (kNearestPoints, currentPoint, lastPoint) {
        if (!lastPoint || lastPoint === currentPoint) {
          lastPoint = [this.maxX.x + 10, this.coords[currentPoint + 1]];
        } else {
          lastPoint = [this.coords[lastPoint], this.coords[lastPoint + 1]];
        }
        this.ray = {
          x0: lastPoint[0],
          y0: lastPoint[1],
          x1: this.coords[currentPoint],
          y1: this.coords[currentPoint + 1]
        };
        const currentPointArr = [this.coords[currentPoint], this.coords[currentPoint + 1]];
        // cant use max or min value for first point, the reference point needs to be the last point in the hull in order to get the angle sorting right
        const rv = sortHeap(this.coords, kNearestPoints, 'polar', lastPoint, currentPointArr).slice();
        // if two points are on the same line eq as current point, currently the further one is considered a 'closer angle', perform swap of these coords below

        let lastSlope;
        let lastDist;
        // if two points relative to each other are in line
        // Issue here when 3 points line up and one is segment from origin
        for (let k = 0; k < rv.length; k++) {
          let lastPoint = [this.coords[rv[k - 1]], this.coords[rv[k - 1] + 1]];
          if (k === 0) {
            lastPoint = currentPointArr;
          }
          const newPoint = [this.coords[rv[k]], this.coords[rv[k] + 1]];
          const newSlope = slope(lastPoint, newPoint);
          const newDist = euclid(currentPointArr, newPoint);

          if (!isNaN(lastSlope) && !isNaN(lastDist) && (Math.abs(newSlope) === Math.abs(lastSlope) || (newSlope === Infinity && lastSlope === -Infinity)) && newDist < lastDist) {
            // flipflop the two points in array order if the slopes are the same
            // sort by euclid instead of straight swap
            swap$1(rv, k, k - 1);
            lastDist = euclid(currentPointArr, [this.coords[rv[k]], this.coords[rv[k] + 1]]);
          } else {
            lastDist = newDist;
          }

          lastSlope = slope(lastPoint, newPoint);
        }

        return rv
      }

      nearestPoints (index, cP, kk) {
        const currentPoint = [this.coords[cP], this.coords[cP + 1]];
        index = sortHeap(this.coords.slice(), index.slice(), 'euclid', currentPoint);
        const rv = [];
        let lastSlope;
        let lastDist;
        kk = Math.min(kk, index.length - 1);
        let i = 0;
        let c = 0;
        if (index.length === 1) {
          return index
        }
        while (c < kk) {
          const newSlope = slope(currentPoint, [this.coords[index[i]], this.coords[index[i] + 1]]);
          const newDist = euclid(currentPoint, [this.coords[index[i]], this.coords[index[i] + 1]]);
          if (newDist <= this.dist) {
            if (newDist === lastDist) {
              rv.push(index[i]);
            } else if (newSlope === lastSlope) {
              const temp = rv[rv.length - 1];
              rv[rv.length - 1] = index[i];
              rv.push(temp);
              c++;
            } else if (c === 0 || (newSlope !== lastSlope)) {
            // } else if ( !lastSlope || newSlope !== lastSlope) {
              rv.push(index[i]);
              c++;
            }
          }
          i++;
          if (i > index.length - 1) {
            return rv
          }
          lastSlope = newSlope;
          lastDist = newDist;
        }

        return rv
      }

      sortHeapAndClean (arr, ind, criteria, minPoint, centerPoint) {
        ind = sortHeap(arr.slice(), ind.slice(), criteria, minPoint, centerPoint);
        ind = this.clean(ind);
        return ind
      }

      clean (index) {
        // there has to be a better way to do this
        // On^2  urrrgh

        let count = 0;
        const duplicates = [];
        for (const item of index) {
          for (let i = 0; i < index.length; i++) {
            if (this.coords[index[i]] === this.coords[item] &&
              this.coords[index[i] + 1] === this.coords[item + 1] &&
              count !== i) {
              let pass = true;
              for (const t of duplicates) {
                if (this.coords[index[i]] === this.coords[t[0]] && this.coords[index[i] + 1] === this.coords[t[0] + 1]) {
                  pass = false;
                  t[1]++;
                }
              }
              if (pass) { duplicates.push([item, 0]); }
              break
            }
          }
          count++;
        }
        const newIndex = [];
        for (let i = 0; i < index.length; i++) {
          let pass = true;
          for (const item of duplicates) {
            if (this.coords[index[i]] === this.coords[item[0]] &&
              this.coords[index[i] + 1] === this.coords[item[0] + 1] &&
              item[1] > 0) {
              item[1]--;
              pass = false;
            }
          }
          if (pass) { newIndex.push(index[i]); }
        }
        // console.log('items removed: ' + (itRem - newIndex.length))
        return newIndex
      }

      calcCenter () {
        const p = { x: 0, y: 0 };

        for (let i = 0; i < this.coords.length; i += 2) {
          p.x += this.coords[i];
          p.y += this.coords[i + 1];
        }
        p.x /= (this.coords.length / 2);
        p.y /= (this.coords.length / 2);
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
        const x1 = point[0] + this.maxD * Math.cos(angle * Math.PI / 180);
        const y1 = point[1] + this.maxD * Math.sin(angle * Math.PI / 180);
        const p = {
          x0: point[0], y0: point[1], x1: x1, y1: y1
        };
        this.ray = p;
        // lets use non-zero winding number rule
        let windingNum = 0;
        let last = { x: Infinity, y: Infinity };

        for (let i = 0; i < index.length; i++) {
          const l = {
            x0: this.coords[index[i]],
            y0: this.coords[index[i] + 1],
            x1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1]],
            y1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1] + 1]
          };
          const inters = intersect(p, l, true);
          if (isFinite(inters.x)) {
            // fails on corner intersection sometimes, so if a corner is intersect, rotate the ray being tested by 1 deg and start over
            const testCond = Math.round(inters.x * 1000000) === last.x && Math.round(inters.y * 1000000) === last.y;
            if (l.y1 - l.y0 > 0 && !testCond) {
              windingNum++;
            } else if (l.y1 - l.y0 < 0 && !testCond) {
              windingNum--;
            } else if (testCond) {
              if (angle < 360) {
                return this.pointInOrOut(point, index, ++angle)
              }
            }
            last = { x: Math.round(inters.x * 1000000), y: Math.round(inters.y * 1000000) };
          }
        }
        return Math.abs(windingNum) !== 0
      }

      printPoints (xIndex) {
        const p = [];
        for (const i of xIndex) {
          p.push(this.coords[i], this.coords[i + 1]);
        }
        console.log(p);
      }

      get hullCoords () {
        return this.subset(this.hull)
      }

      subset (indices) {
        const rv = [];
        for (const i of indices) {
          rv.push(this.coords[i], this.coords[i + 1]);
        }
        return rv
      }

      get coords2D () {
        const newArr = [];
        const arr = this.sortedCoords;
        while (arr.length) newArr.push(arr.splice(0, 2));
        return newArr
      }

      get sortedCoords () {
        const newArr = [];
        for (const i of this.hull) {
          newArr.push(this.coords[i], this.coords[i + 1]);
        }
        return newArr
      }
    }

    class BoundaryExtra extends Boundary {
      constructor (arr, k = 3, maxDist = Infinity) {
        super(arr, k, maxDist);
        //    this.cPoints = []
        this.origCoordsLen = arr.length;
        this.intersectingLineSegs = [];
        this.maxDist = maxDist;
      }

      /**
       * addPoints
       * use final k value from concave boundary for point search in order to
       * @param {Array} parentArr Array of coordinate cloud used to interpolate boundary points
       * @param {Object} constructed delaunator object, need this for edges and triangles
       * @param {Integer} dist Max distance to point to trigger interpolation, only one of two points in line segment has to meet this criteria
       */
      addPoints (parentArr, delaunator, dist) {
        this.k = 3;
        const edges = getEdges(delaunator);
        // get all intersecting lines to the hull line seg
        for (let p = 0; p < this.hull.length - 1; p++) {
          const h = this.subset([this.hull[p], this.hull[p + 1]]);
          const seg = { x0: h[0], y0: h[1], x1: h[2], y1: h[3] };
          const temp = this.getIntersectingLines(seg, edges, parentArr, dist).reverse();
          const ind = sortHeap(temp.map((m) => [m[m.length - 1].x, m[m.length - 1].y]).flat(), [...Array(temp.length).keys()].map((i) => i * 2), 'euclid', [seg.x0, seg.y0]);
          this.intersectingLineSegs.push(this.hull[p]);
          const c = this.coords.length;
          this.coords = this.coords.concat(temp.map((m) => [m[m.length - 1].x, m[m.length - 1].y]).flat());
          this.intersectingLineSegs = this.intersectingLineSegs.concat(ind.map((m) => m + c));
        }
        this.intersectingLineSegs.push(this.intersectingLineSegs[0]);

        this.hull = this.intersectingLineSegs;
        this.clean(this.hull);
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
        const pntAndItsArr = [];
        // iterate over point pairs
        for (let i = 0; i < indexArr.length; i += 2) {
          // first check if distance from either point to perpendicular of line seg is less than dist
          const point = { x: coords[indexArr[i]], y: coords[indexArr[i] + 1] };
          const d = distLineAndPoint(lineSeg, point);
          if (d <= dist) {
            if (i + 1 > indexArr.length - 1) {
              i = -1;
            }
            const coordSeg = { x0: point.x, y0: point.y, x1: coords[indexArr[i + 1]], y1: coords[indexArr[i + 1] + 1] };
            const its = intersect(lineSeg, coordSeg, false);
            if (isFinite(its.x)) {
              // console.log(d, its, lineSeg, coordSeg, i)
              // this.cPoints.push(indexArr[i], indexArr[i + 1])
              if (opt) {
                opt([indexArr[i], indexArr[i + 1], its]);
              }
              pntAndItsArr.push([indexArr[i], indexArr[i + 1], its]);
            }
          }
        }
        return pntAndItsArr
      }

      get k () {
        return super.k
      }

      set k (v) {
        super.k = v;
      }
    }

    /**
     * ConstrainoDelaunato
     *
     * @class
     * @classdesc ConstrainoDelaunato
     */
    class ConstrainoDelaunato {
      /**
       * creates a delaunator object for the larger coord point cloud, and any smalle concave boundaries and delaunator objects for holes/boundaries supplied
       * @constructor
       *
       * @param {Array} coords Coordinate cloud, can be 2D or 1D, prefer 1D of type [x0, y0, x1, y1, ... xN, yN]
       * @param {Integer} k lower bound for point selection in k grouping - minimum possible value is 3 - you have to make a polygon
       * @param {Integer} dist - distance for adding points along boundary, distance between line segment perpindicular to either point of triangle segment; used for point interpolation along a boundary
       * @param {Integer} distSelectionLimit - distance to limit selection of candidate points for concave boundary creation - useful if there is a whole on edge of points that is not being acknowledged by algorithm due to uniform point spacing or something like that; used during concave boundary creation
       * @param {Array} ...boundaries Point clouds of holes in coords, stored in array boundary for concave boundaries and boundedDelaunator for created delaunator objects
       */
      constructor (coords, k, calcDist, ...boundaries) {
        // k is the k-nearest neighbor selection
        // if coords are 2D
        if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
          coords = coords.flat();
        } else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
          return
        }

        this._delaunator = new Delaunator(coords);
        const edges = getEdges(this._delaunator);
        // calculate a average dist as the max selection distance
        let count = 0;
        let sum = 0;
        for (let i = 0; i < edges.length; i += 2) {
          if (i === edges.length - 1) {
            sum += euclid([coords[edges[i]], coords[edges[i] + 1]], [coords[edges[0]], coords[edges[0] + 1]]);
          } else {
            sum += euclid([coords[edges[i]], coords[edges[i] + 1]], [coords[edges[i + 1]], coords[edges[i + 1] + 1]]);
          }
          count += 1;
        }

        let distSelectionLimit = Infinity;
        if (calcDist) {
          distSelectionLimit = Math.ceil(sum / count);
        }

        const dist = Math.ceil(sum / count) * 2;
        this._boundaries = [];
        this.boundedDelaunators = [];

        for (let boundary of boundaries) {
          if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
            boundary = boundary.flat();
          }
          this._boundaries.push(new BoundaryExtra(boundary, k));
          this._boundaries[this._boundaries.length - 1].addPoints(coords, this._delaunator, dist);
          this.boundedDelaunators.push(this.setTrianglesInsideBound(this._boundaries[this._boundaries.length - 1]));
        }
        this._boundaries.push(new BoundaryExtra(coords, k, distSelectionLimit));
        this._boundaries[this._boundaries.length - 1].addPoints(coords, this._delaunator, dist);
        this.boundedDelaunators.push(this.setTrianglesInsideBound(this._boundaries[this._boundaries.length - 1]));
      }

      /**
       * setTrianglesInsideBound
       *
       * Function used to clip coords to inside of boundary or hole
       *
       * @param {BoundaryExtra} boundary boundary extra object
       */
      setTrianglesInsideBound (boundary) {
        let coords = [];
        const index = [...this.delaunator.coords.keys()].filter((i) => i % 2 === 0);
        const maxX = maximumPointX(this.delaunator.coords, index);
        for (const e of index) {
          const point = { x: this.delaunator.coords[e], y: this.delaunator.coords[e + 1] };
            coords.push(point.x, point.y);
          // if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
          //   outIndex.push(i)
          //   coords.push(point.x, point.y)
          //   i += 2
          // }
        }

        console.log(coords.length, index, coords);
        coords = coords.concat(boundary.subset(boundary.hull));
        const rv = new Delaunator(coords);
        const t = [];
        for (let e = 0; e < rv.triangles.length / 3; e++) {
          const edgeIndex = e * 3;

          let xCoord = 0;
          let yCoord = 0;
          for (let r = 0; r < 3; r++) {
            xCoord += rv.coords[2 * rv.triangles[r + edgeIndex]];
            yCoord += rv.coords[2 * rv.triangles[r + edgeIndex] + 1];
          }
          const point = { x: xCoord / 3, y: yCoord / 3 };
          if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
            t.push(rv.triangles[edgeIndex], rv.triangles[edgeIndex + 1], rv.triangles[edgeIndex + 2]);
          }
        }
        rv.triangles = new rv.triangles.constructor(t);

        return rv
      }

      /**
       * update
       *
       * @param {Array} point x and y coord of point to add the delaunator object
       */
      update (point) {
        const c = this.coords;
        for (const p of point.flat()) {
          c.push(p);
        }
        this._delaunator = new Delaunator(c);
      }

      /**
       * coords2D
       *
       * @returns {Array} 2D coordinate array
       */
      get coords2D () {
        const c2D = [];
        const c1D = this.coords;
        for (let i = 0; i < c1D.length; i += 2) {
          c2D.push([c1D[i], c1D[i + 1]]);
        }
        return c2D
      }

      /**
       * coords
       *
       * @returns {Array} 1D coordinate array
       */
      get coords () {
        return this.delaunator.coords
      }

      /**
       * triangles
       *
       * @returns {Array} Index array of delaunator triangles
       */
      get triangles () {
        return this.delaunator.triangles
      }

      /**
       * hull
       *
       * @returns {Array} Array of hull indices
       */
      get hull () {
        return this.delaunator.hull
      }

      /**
       * delaunator
       *
       * @returns {Object} Return delaunator object for parent points
       */
      get delaunator () {
        return this._delaunator
      }

      /**
       * boundaries
       *
       * @returns {Array} concave boundary index array of parent points and hole points includes parent points as final array item
       */
      get boundaries () {
        return this._boundaries
      }

      /**
       * holes
       *
       * @returns {Array} Delaunator object array for all hole/boundary points supplied, includes parent points as final array item
       */
      get holes () {
        return this.boundedDelaunators
      }
    }

    return ConstrainoDelaunato;

})));
