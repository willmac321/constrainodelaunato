(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global = global || self, global.ConstrainoDelaunato = factory());
}(this, (function () { 'use strict';

    var EPSILON = Math.pow(2, -52);
    var EDGE_STACK = new Uint32Array(512);

    var Delaunator = function Delaunator(coords) {
        var n = coords.length >> 1;
        if (n > 0 && typeof coords[0] !== 'number') { throw new Error('Expected coords to contain numbers.'); }

        this.coords = coords;

        // arrays that will store the triangulation graph
        var maxTriangles = Math.max(2 * n - 5, 0);
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
    };

    Delaunator.from = function from (points, getX, getY) {
            if ( getX === void 0 ) getX = defaultGetX;
            if ( getY === void 0 ) getY = defaultGetY;

        var n = points.length;
        var coords = new Float64Array(n * 2);

        for (var i = 0; i < n; i++) {
            var p = points[i];
            coords[2 * i] = getX(p);
            coords[2 * i + 1] = getY(p);
        }

        return new Delaunator(coords);
    };

    Delaunator.prototype.update = function update () {
        var ref =  this;
            var coords = ref.coords;
            var hullPrev = ref._hullPrev;
            var hullNext = ref._hullNext;
            var hullTri = ref._hullTri;
            var hullHash = ref._hullHash;
        var n = coords.length >> 1;

        // populate an array of point indices; calculate input data bbox
        var minX = Infinity;
        var minY = Infinity;
        var maxX = -Infinity;
        var maxY = -Infinity;

        for (var i = 0; i < n; i++) {
            var x = coords[2 * i];
            var y = coords[2 * i + 1];
            if (x < minX) { minX = x; }
            if (y < minY) { minY = y; }
            if (x > maxX) { maxX = x; }
            if (y > maxY) { maxY = y; }
            this._ids[i] = i;
        }
        var cx = (minX + maxX) / 2;
        var cy = (minY + maxY) / 2;

        var minDist = Infinity;
        var i0, i1, i2;

        // pick a seed point close to the center
        for (var i$1 = 0; i$1 < n; i$1++) {
            var d = dist(cx, cy, coords[2 * i$1], coords[2 * i$1 + 1]);
            if (d < minDist) {
                i0 = i$1;
                minDist = d;
            }
        }
        var i0x = coords[2 * i0];
        var i0y = coords[2 * i0 + 1];

        minDist = Infinity;

        // find the point closest to the seed
        for (var i$2 = 0; i$2 < n; i$2++) {
            if (i$2 === i0) { continue; }
            var d$1 = dist(i0x, i0y, coords[2 * i$2], coords[2 * i$2 + 1]);
            if (d$1 < minDist && d$1 > 0) {
                i1 = i$2;
                minDist = d$1;
            }
        }
        var i1x = coords[2 * i1];
        var i1y = coords[2 * i1 + 1];

        var minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (var i$3 = 0; i$3 < n; i$3++) {
            if (i$3 === i0 || i$3 === i1) { continue; }
            var r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i$3], coords[2 * i$3 + 1]);
            if (r < minRadius) {
                i2 = i$3;
                minRadius = r;
            }
        }
        var i2x = coords[2 * i2];
        var i2y = coords[2 * i2 + 1];

        if (minRadius === Infinity) {
            // order collinear points by dx (or dy if all x are identical)
            // and return the list as a hull
            for (var i$4 = 0; i$4 < n; i$4++) {
                this._dists[i$4] = (coords[2 * i$4] - coords[0]) || (coords[2 * i$4 + 1] - coords[1]);
            }
            quicksort(this._ids, this._dists, 0, n - 1);
            var hull = new Uint32Array(n);
            var j = 0;
            for (var i$5 = 0, d0 = -Infinity; i$5 < n; i$5++) {
                var id = this._ids[i$5];
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
            var i$6 = i1;
            var x$1 = i1x;
            var y$1 = i1y;
            i1 = i2;
            i1x = i2x;
            i1y = i2y;
            i2 = i$6;
            i2x = x$1;
            i2y = y$1;
        }

        var center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
        this._cx = center.x;
        this._cy = center.y;

        for (var i$7 = 0; i$7 < n; i$7++) {
            this._dists[i$7] = dist(coords[2 * i$7], coords[2 * i$7 + 1], center.x, center.y);
        }

        // sort the points by distance from the seed triangle circumcenter
        quicksort(this._ids, this._dists, 0, n - 1);

        // set up the seed triangle as the starting hull
        this._hullStart = i0;
        var hullSize = 3;

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

        for (var k = 0, xp = (void 0), yp = (void 0); k < this._ids.length; k++) {
            var i$8 = this._ids[k];
            var x$2 = coords[2 * i$8];
            var y$2 = coords[2 * i$8 + 1];

            // skip near-duplicate points
            if (k > 0 && Math.abs(x$2 - xp) <= EPSILON && Math.abs(y$2 - yp) <= EPSILON) { continue; }
            xp = x$2;
            yp = y$2;

            // skip seed triangle points
            if (i$8 === i0 || i$8 === i1 || i$8 === i2) { continue; }

            // find a visible edge on the convex hull using edge hash
            var start = 0;
            for (var j$1 = 0, key = this._hashKey(x$2, y$2); j$1 < this._hashSize; j$1++) {
                start = hullHash[(key + j$1) % this._hashSize];
                if (start !== -1 && start !== hullNext[start]) { break; }
            }

            start = hullPrev[start];
            var e = start, q = (void 0);
            while (q = hullNext[e], !orient(x$2, y$2, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
                e = q;
                if (e === start) {
                    e = -1;
                    break;
                }
            }
            if (e === -1) { continue; } // likely a near-duplicate point; skip it

            // add the first triangle from the point
            var t = this._addTriangle(e, i$8, hullNext[e], -1, -1, hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i$8] = this._legalize(t + 2);
            hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize++;

            // walk forward through the hull, adding more triangles and flipping recursively
            var n$1 = hullNext[e];
            while (q = hullNext[n$1], orient(x$2, y$2, coords[2 * n$1], coords[2 * n$1 + 1], coords[2 * q], coords[2 * q + 1])) {
                t = this._addTriangle(n$1, i$8, q, hullTri[i$8], -1, hullTri[n$1]);
                hullTri[i$8] = this._legalize(t + 2);
                hullNext[n$1] = n$1; // mark as removed
                hullSize--;
                n$1 = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e === start) {
                while (q = hullPrev[e], orient(x$2, y$2, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
                    t = this._addTriangle(q, i$8, e, -1, hullTri[e], hullTri[q]);
                    this._legalize(t + 2);
                    hullTri[q] = t;
                    hullNext[e] = e; // mark as removed
                    hullSize--;
                    e = q;
                }
            }

            // update the hull indices
            this._hullStart = hullPrev[i$8] = e;
            hullNext[e] = hullPrev[n$1] = i$8;
            hullNext[i$8] = n$1;

            // save the two new edges in the hash table
            hullHash[this._hashKey(x$2, y$2)] = i$8;
            hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
        }

        this.hull = new Uint32Array(hullSize);
        for (var i$9 = 0, e$1 = this._hullStart; i$9 < hullSize; i$9++) {
            this.hull[i$9] = e$1;
            e$1 = hullNext[e$1];
        }

        // trim typed triangle mesh arrays
        this.triangles = this._triangles.subarray(0, this.trianglesLen);
        this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
    };

    Delaunator.prototype._hashKey = function _hashKey (x, y) {
        return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
    };

    Delaunator.prototype._legalize = function _legalize (a) {
        var ref = this;
            var triangles = ref._triangles;
            var halfedges = ref._halfedges;
            var coords = ref.coords;

        var i = 0;
        var ar = 0;

        // recursion eliminated with a fixed-size stack
        while (true) {
            var b = halfedges[a];

            /* if the pair of triangles doesn't satisfy the Delaunay condition
             * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
             * then do the same check/flip recursively for the new pair of triangles
             *
             *       pl                pl
             *      /||\              /  \
             *   al/ || \bl        al/\a
             *    /  ||  \          /  \
             *   /  a||b  \flip/___ar___\
             * p0\   ||   /p1   =>   p0\---bl---/p1
             *    \  ||  /          \  /
             *   ar\ || /br         b\/br
             *      \||/              \  /
             *       pr                pr
             */
            var a0 = a - a % 3;
            ar = a0 + (a + 2) % 3;

            if (b === -1) { // convex hull edge
                if (i === 0) { break; }
                a = EDGE_STACK[--i];
                continue;
            }

            var b0 = b - b % 3;
            var al = a0 + (a + 1) % 3;
            var bl = b0 + (b + 2) % 3;

            var p0 = triangles[ar];
            var pr = triangles[a];
            var pl = triangles[al];
            var p1 = triangles[bl];

            var illegal = inCircle(
                coords[2 * p0], coords[2 * p0 + 1],
                coords[2 * pr], coords[2 * pr + 1],
                coords[2 * pl], coords[2 * pl + 1],
                coords[2 * p1], coords[2 * p1 + 1]);

            if (illegal) {
                triangles[a] = p1;
                triangles[b] = p0;

                var hbl = halfedges[bl];

                // edge swapped on the other side of the hull (rare); fix the halfedge reference
                if (hbl === -1) {
                    var e = this._hullStart;
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

                var br = b0 + (b + 1) % 3;

                // don't worry about hitting the cap: it can only happen on extremely degenerate input
                if (i < EDGE_STACK.length) {
                    EDGE_STACK[i++] = br;
                }
            } else {
                if (i === 0) { break; }
                a = EDGE_STACK[--i];
            }
        }

        return ar;
    };

    Delaunator.prototype._link = function _link (a, b) {
        this._halfedges[a] = b;
        if (b !== -1) { this._halfedges[b] = a; }
    };

    // add a new triangle given vertex indices and adjacent half-edge ids
    Delaunator.prototype._addTriangle = function _addTriangle (i0, i1, i2, a, b, c) {
        var t = this.trianglesLen;

        this._triangles[t] = i0;
        this._triangles[t + 1] = i1;
        this._triangles[t + 2] = i2;

        this._link(t, a);
        this._link(t + 1, b);
        this._link(t + 2, c);

        this.trianglesLen += 3;

        return t;
    };

    // monotonically increases with real angle, but doesn't need expensive trigonometry
    function pseudoAngle(dx, dy) {
        var p = dx / (Math.abs(dx) + Math.abs(dy));
        return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
    }

    function dist(ax, ay, bx, by) {
        var dx = ax - bx;
        var dy = ay - by;
        return dx * dx + dy * dy;
    }

    // return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
    function orientIfSure(px, py, rx, ry, qx, qy) {
        var l = (ry - py) * (qx - px);
        var r = (rx - px) * (qy - py);
        return Math.abs(l - r) >= 3.3306690738754716e-16 * Math.abs(l + r) ? l - r : 0;
    }

    // a more robust orientation test that's stable in a given triangle (to fix robustness issues)
    function orient(rx, ry, qx, qy, px, py) {
        var sign = orientIfSure(px, py, rx, ry, qx, qy) ||
        orientIfSure(rx, ry, qx, qy, px, py) ||
        orientIfSure(qx, qy, px, py, rx, ry);
        return sign < 0;
    }

    function inCircle(ax, ay, bx, by, cx, cy, px, py) {
        var dx = ax - px;
        var dy = ay - py;
        var ex = bx - px;
        var ey = by - py;
        var fx = cx - px;
        var fy = cy - py;

        var ap = dx * dx + dy * dy;
        var bp = ex * ex + ey * ey;
        var cp = fx * fx + fy * fy;

        return dx * (ey * cp - bp * fy) -
               dy * (ex * cp - bp * fx) +
               ap * (ex * fy - ey * fx) < 0;
    }

    function circumradius(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = 0.5 / (dx * ey - dy * ex);

        var x = (ey * bl - dy * cl) * d;
        var y = (dx * cl - ex * bl) * d;

        return x * x + y * y;
    }

    function circumcenter(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = 0.5 / (dx * ey - dy * ex);

        var x = ax + (ey * bl - dy * cl) * d;
        var y = ay + (dx * cl - ex * bl) * d;

        return {x: x, y: y};
    }

    function quicksort(ids, dists, left, right) {
        if (right - left <= 20) {
            for (var i = left + 1; i <= right; i++) {
                var temp = ids[i];
                var tempDist = dists[temp];
                var j = i - 1;
                while (j >= left && dists[ids[j]] > tempDist) { ids[j + 1] = ids[j--]; }
                ids[j + 1] = temp;
            }
        } else {
            var median = (left + right) >> 1;
            var i$1 = left + 1;
            var j$1 = right;
            swap(ids, median, i$1);
            if (dists[ids[left]] > dists[ids[right]]) { swap(ids, left, right); }
            if (dists[ids[i$1]] > dists[ids[right]]) { swap(ids, i$1, right); }
            if (dists[ids[left]] > dists[ids[i$1]]) { swap(ids, left, i$1); }

            var temp$1 = ids[i$1];
            var tempDist$1 = dists[temp$1];
            while (true) {
                do { i$1++; } while (dists[ids[i$1]] < tempDist$1);
                do { j$1--; } while (dists[ids[j$1]] > tempDist$1);
                if (j$1 < i$1) { break; }
                swap(ids, i$1, j$1);
            }
            ids[left + 1] = ids[j$1];
            ids[j$1] = temp$1;

            if (right - i$1 + 1 >= j$1 - left) {
                quicksort(ids, dists, i$1, right);
                quicksort(ids, dists, left, j$1 - 1);
            } else {
                quicksort(ids, dists, left, j$1 - 1);
                quicksort(ids, dists, i$1, right);
            }
        }
    }

    function swap(arr, i, j) {
        var tmp = arr[i];
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
      var rv = [];
      for (var e = 0; e < delaunay.triangles.length; e++) {
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
    function intersect (p, l, checkEndpoints) {
      if ( checkEndpoints === void 0 ) checkEndpoints = false;

      // compare two line segments to see if they intersect
      var den = ((l.y1 - l.y0) * (p.x1 - p.x0)) - ((l.x1 - l.x0) * (p.y1 - p.y0));
      if (den === 0) {
        return { x: Infinity, y: Infinity }
      }

      var a = p.y0 - l.y0;
      var b = p.x0 - l.x0;

      var num1 = ((l.x1 - l.x0) * a) - ((l.y1 - l.y0) * b);
      var num2 = ((p.x1 - p.x0) * a) - ((p.y1 - p.y0) * b);

      a = num1 / den;
      b = num2 / den;

      var rv = {
        x: p.x0 + (a * (p.x1 - p.x0)),
        y: p.y0 + (a * (p.y1 - p.y0))
      };

      // if (p.y1 === rv.y) {
      // console.log(a, b);
      // }
      //
      var t = compareIntersect(a, b);
      if (checkEndpoints) {
        t = compareIntersectEndpoints(a, b);
      }

      if (t.a && t.b) {
        return rv
      }
      return { x: Infinity, y: Infinity }
    }

    function dotProduct (a, b) {
      var p = { x: a[0], y: a[1] };
      var o = { x: b[0], y: b[1] };
      return p.x * o.x + p.y * o.y
    }

    function slope (a, b) {
      var p = { x: a[0], y: a[1] };
      var o = { x: b[0], y: b[1] };
      return (p.y - o.y) / (p.x - o.x)
    }

    function compareIntersect (a, b) {
      var t = { a: false, b: false };
      if (a > 0 && a < 1) {
        t.a = true;
      }
      if (b > 0 && b < 1) {
        t.b = true;
      }
      return t
    }

    function compareIntersectEndpoints (a, b) {
      var t = { a: false, b: false };
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
      var o = { x: a[0] - c.x, y: a[1] - c.y };
      var p = { x: b[0] - c.x, y: b[1] - c.y };
      var magO = Math.sqrt(Math.pow(o.x, 2) + Math.pow(o.y, 2));
      var magP = Math.sqrt(Math.pow(p.x, 2) + Math.pow(p.y, 2));
      var theta = Math.acos((p.x * o.x + p.y * o.y) / (magO * magP)) * 180 / Math.PI;
      var det = (o.x * p.y) - (p.x * o.y);
      theta = isNaN(theta) ? 0 : theta;
      theta = det > 0 ? 360 - theta : theta;
      return theta
    }

    function manhattenDist (a, b) {
      var p = { x: a[0], y: a[1] };
      var o = { x: b[0], y: b[1] };
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
      var p = { x: a[0], y: a[1] };
      var o = { x: b[0], y: b[1] };
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
      var num = Math.abs((l.y1 - l.y0) * p.x - (l.x1 - l.x0) * p.y + l.x1 * l.y0 - l.y1 * l.x0);
      var dem = Math.sqrt(Math.pow(l.y1 - l.y0, 2) + Math.pow(l.x1 - l.x0, 2));
      return num / dem
    }

    // heap sort 2d array by angle
    function heapSort (minpoint, index, a, count, p, center) {
      var func;

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
        func = function (a, b) { return a - b; };
      } else {
        func = dotProduct;
      }

      heapify(minpoint, a, index, count, func, p);

      var end = count - 1;
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
      var par = function (i) { return Math.floor((i - 1) / 2); };
      var start = par(count - 1);
      while (start >= 0) {
        siftDown(point, a, index, start, count - 1, func, p);
        start--;
      }
    }

    function siftDown (point, a, index, start, end, func, p) {
      var root = start;
      var left = function (i) { return 2 * i + 1; };

      while (left(root) <= end) {
        var child = left(root);
        var s = root;
        var dot = function (i) {
          var t = func([a[index[i]], a[index[i] + 1]], point);
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
      var t = a[i];
      a[i] = a[j];
      a[j] = t;
    }

    function maximumPointX (newArr, index) {
      var ind = 0;
      var minY = -Infinity;
      var minX = -Infinity;
      if (index) {
        for (var i = 0, list = index.entries(); i < list.length; i += 1) {
          var ref = list[i];
          var k = ref[0];
          var p = ref[1];

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
        for (var p$1 = 0; p$1 > newArr.length; p$1++) {
          if (newArr[p$1] > minX) {
            minX = newArr[p$1];
            ind = p$1;
          }
        }
      }
      return { x: minX, y: minY, i: ind }
    }

    function minimumPointY (newArr, index) {
      var ind = 0;
      var minY = Infinity;
      var minX = Infinity;
      if (index) {
        for (var i = 0, list = index.entries(); i < list.length; i += 1) {
          var ref = list[i];
          var k = ref[0];
          var p = ref[1];

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
        for (var p$1 = 0; p$1 < newArr.length; p$1++) {
          if (newArr[p$1] < minX) {
            minX = newArr[p$1];
            ind = p$1;
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
      var ind = 0;
      var minY = Infinity;
      var minX = Infinity;
      if (index) {
        for (var i = 0, list = index.entries(); i < list.length; i += 1) {
          var ref = list[i];
          var k = ref[0];
          var p = ref[1];

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
        for (var p$1 = 0; p$1 < newArr.length; p$1++) {
          if (newArr[p$1] < minX) {
            minX = newArr[p$1];
            ind = p$1;
          }
        }
      }
      return { x: minX, y: minY, i: ind }
    }

    var counter = 0;

    var Boundary = function Boundary (arr, k) {
      if ( k === void 0 ) k = 3;

      this.k = k;
      this.coords = arr.slice();
      this.index = [].concat( this.coords.keys() ).filter(function (i) { return i % 2 === 0; });
      this.index = this.clean(this.index);
      this.center = this.calcCenter();
      this.minY = minimumPointY(this.coords, this.index);
      this.minX = minimumPointX(this.coords, this.index);
      this.maxX = maximumPointX(this.coords, this.index);

      this.cPoints = [];

      this.ray = null;
      this.hull = this.findConcaveHull(k);
    };

    var prototypeAccessors = { hullCoords: { configurable: true },coords2D: { configurable: true },sortedCoords: { configurable: true } };

    /**
     * findConcaveHull
     *
     * @param {Integer} k Starting point cloud count for points being used for selection
     * @returns {Array} Array of indices of X values for the sorted concave hull
     */
    Boundary.prototype.findConcaveHull = function findConcaveHull (k) {
      // alt index is sorted to minX value
      var index = this.sortHeapAndClean(this.coords, this.index, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y]);
      var hull = this.concave(index, k);
      // hull = this.sortHeapAndClean(this.coords, hull, 'polar', [this.minX.x, this.minY.y], [this.center.x, this.center.y])
      // hull.push(hull[0])
      return hull
    };

    Boundary.prototype.concave = function concave (index, k) {
      // k nearest neighbor babbbbyyyy
      // https://pdfs.semanticscholar.org/2397/17005c3ebd5d6a42fc833daf97a0edee1ce4.pdf
      // double check arr is sorted and clean
      // also sort it so all points are in order from some min pointon the xy plane
      var stopVal = Infinity; // 76 // Infinity // and beyond
      var oldIndex = index.slice();
      // console.log('new k', k)
      if (index.length < 3) {
        console.log('len less than 3');
        return null
      } else if (k > index.length - 1) {
        console.log(counter);
        console.log('k is too big');
        return null
      } else if (index.length === 3) {
        console.log('len 3');
        return index
      }

      var kk = Math.min(Math.max(k, 3), index.length - 1);
      // i is a pointer to the relative index not a loc in this.coords
      // so, index of that index gives a this.coords pointer
      var firstPointIndex = minimumPointY(this.coords, index).i;
      var firstPoint = { i: firstPointIndex, coord: index[firstPointIndex] };
      var currentPoint = firstPoint.coord;
      var hull = [firstPoint.coord];
      // why is step init to 2?
      // Because the paper was written in Matlab....
      var step = 1;
      // each index value can only be used once so this is ok
      index.splice(firstPoint.i, 1);
      while ((currentPoint !== firstPoint.coord || step === 1) && (index.length > 0)) {
        counter++;
        if (step === 4) {
          index.push(firstPoint.coord);
        }
        // find nearest neighbors
        var kNearestPoints = this.nearestPoints(index, currentPoint, kk);
        // descending order 'right-hand' turn x and y min are top left on js canvas in webpage
        var cPoints = this.sortByAngle(kNearestPoints, currentPoint, hull[hull.length - 2]);
        // if (cPoints.indexOf(firstPoint.coord) > -1) {
        // console.log(cPoints)
        // }
        var its = true;
        var i = -1;
        while (its && i < cPoints.length - 1) {
          // This is so that when the first point is added to the end of the hull, it doesn't get used to check for intersections
          var lastPoint = 0;
          if (cPoints[i] === firstPoint.coord) {
            // console.log('back to first', firstPoint)
            lastPoint = 1;
          }
          var j = 1;
          its = false;
          while (!its && j < hull.length - lastPoint) {
            var l = {
              x0: this.coords[hull[step - 1]],
              y0: this.coords[hull[step - 1] + 1],
              x1: this.coords[cPoints[i + 1]],
              y1: this.coords[cPoints[i + 1] + 1]
            };
            var p = {
              x0: this.coords[hull[step - j]],
              y0: this.coords[hull[step - j] + 1],
              x1: this.coords[hull[step - 1 - j]],
              y1: this.coords[hull[step - 1 - j] + 1]
            };
            // the endpoint of one line segment is always intersecting the endpoint of a connected line segment, how to ignore this intersection?
            var ints = intersect(p, l, true);
            var endpointsMatch = (p.x0 === l.x0 && p.y0 === l.y0);
            var isClose = (cPoints[i + 1] === firstPoint.coord) && (p.x1 === l.x1 && p.y1 === l.y1);
            // (p.x0 !== l.x0 && p.y0 !== l.y0) ||
            // if (l.x0 === 221 && l.y0 === 90) {
            // console.log(l, p, ints, isFinite(ints.x), (p.x1 === l.x0 && p.y1 === l.y0))
            // }
            if (isFinite(ints.x) && !endpointsMatch && !isClose) {
              its = true;
            }
            // if (counter > 269) {
            // console.log(l, p, ints, isFinite(ints.x), !endpointsMatch, cPoints[i + 1], firstPoint)
            // console.log(its)
            // }
            j++;
          }
          i++;
        }
        this.cPoints = cPoints.slice();
        // this.cPoints.splice(i, 1)

        if (its) {
          // console.log('intersection found at k ', k, its)
          // if (kk + 1 === 12) {
          // console.log(counter)
          // return hull
          // }
          return this.concave(oldIndex, ++kk)
        }
        currentPoint = cPoints[i];
        hull.push(currentPoint);

        if (counter > stopVal) {
          // console.log('test', [this.coords[currentPoint], this.coords[currentPoint + 1]], this.subset(cPoints))
          return hull // .concat(cPoints)
        }
        index.splice(index.indexOf(currentPoint), 1);
        step++;
      }
      var allInside = true;
      for (var i$2 = 0, list = index; i$2 < list.length; i$2 += 1) {
        var i$1 = list[i$2];

          allInside = this.pointInOrOut(
          [this.coords[i$1], this.coords[i$1 + 1]],
          hull, this.maxX.x + 10);
        if (!allInside) {
          break
        }
      }
      if (!allInside) {
        // console.log('Another time round')
        //if (kk + 1 === 11) {
        //  console.log(counter)
        //  return hull
        //}
        return this.concave(oldIndex, ++kk)
      }
      // console.log('made it out')
      this.k = kk;
      return hull
    };

    Boundary.prototype.sortByAngle = function sortByAngle (kNearestPoints, currentPoint, lastPoint) {
      // const lastPointIndex = lastPoint
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
      var currentPointArr = [this.coords[currentPoint], this.coords[currentPoint + 1]];
      // cant use max or min value for first point, the reference point needs to be the last point in the hull in order to get the angle sorting right
      var rv = sortHeap(this.coords, kNearestPoints, 'polar', lastPoint, currentPointArr).slice();
      // if two points are on the same line eq as current point, currently the further one is considered a 'closer angle', perform swap of these coords below

      var lastSlope;
      var lastDist;
      // if two points relative to each other are in line
      // Issue here when 3 points line up and one is segment from origin
      for (var k = 0; k < rv.length; k++) {
        var lastPoint$1 = [this.coords[rv[k - 1]], this.coords[rv[k - 1] + 1]];
        if (k === 0) {
          lastPoint$1 = currentPointArr;
        }
        var newPoint = [this.coords[rv[k]], this.coords[rv[k] + 1]];
        var newSlope = slope(lastPoint$1, newPoint);
        var newDist = euclid(currentPointArr, newPoint);
        if (lastSlope && lastDist && (Math.abs(newSlope) === Math.abs(lastSlope) || (newSlope === Infinity && lastSlope === -Infinity)) && newDist < lastDist) {
          // flipflop the two points in array order if the slopes are the same
          // sort by euclid instead of straight swap
          swap$1(rv, k, k - 1);
          lastDist = euclid(currentPointArr, [this.coords[rv[k]], this.coords[rv[k] + 1]]);
        } else {
          lastDist = newDist;
        }

        lastSlope = slope(lastPoint$1, newPoint);
      }

      return rv
    };

    Boundary.prototype.nearestPoints = function nearestPoints (index, cP, kk) {
      var currentPoint = [this.coords[cP], this.coords[cP + 1]];
      index = sortHeap(this.coords.slice(), index.slice(), 'euclid', currentPoint);
      var rv = [];
      var lastSlope;
      kk = Math.min(kk, index.length - 1);
      var i = 0;
      var c = 0;
      while (c < kk) {
        var newSlope = slope(currentPoint, [this.coords[index[i]], this.coords[index[i] + 1]]);
        if (!lastSlope || newSlope !== lastSlope) {
          rv.push(index[i]);
          c++;
        }
        i++;
        if (i > index.length - 1) {
          return rv
        }
        lastSlope = newSlope;
      }
      return rv
    };

    Boundary.prototype.sortHeapAndClean = function sortHeapAndClean (arr, ind, criteria, minPoint, centerPoint) {
      ind = sortHeap(arr.slice(), ind.slice(), criteria, minPoint, centerPoint);
      ind = this.clean(ind);
      return ind
    };

    Boundary.prototype.clean = function clean (index) {
      // there has to be a better way to do this
      // On^2urrrgh
      var itRem = index.length;

      var count = 0;
      var duplicates = [];
      for (var i$3 = 0, list$1 = index; i$3 < list$1.length; i$3 += 1) {
        var item = list$1[i$3];

          for (var i = 0; i < index.length; i++) {
          if (this.coords[index[i]] === this.coords[item] &&
            this.coords[index[i] + 1] === this.coords[item + 1] &&
            count !== i) {
            var pass = true;
            for (var i$2 = 0, list = duplicates; i$2 < list.length; i$2 += 1) {
              var t = list[i$2];

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
      var newIndex = [];
      for (var i$1 = 0; i$1 < index.length; i$1++) {
        var pass$1 = true;
        for (var i$4 = 0, list$2 = duplicates; i$4 < list$2.length; i$4 += 1) {
          var item$1 = list$2[i$4];

            if (this.coords[index[i$1]] === this.coords[item$1[0]] &&
            this.coords[index[i$1] + 1] === this.coords[item$1[0] + 1] &&
            item$1[1] > 0) {
            item$1[1]--;
            pass$1 = false;
          }
        }
        if (pass$1) { newIndex.push(index[i$1]); }
      }
      console.log('items removed: ' + (itRem - newIndex.length));
      return newIndex
    };

    Boundary.prototype.calcCenter = function calcCenter () {
      var p = { x: 0, y: 0 };

      for (var i = 0; i < this.coords.length; i += 2) {
        p.x += this.coords[i];
        p.y += this.coords[i + 1];
      }
      p.x /= (this.coords.length / 2);
      p.y /= (this.coords.length / 2);
      return p
    };

    Boundary.prototype.pointInOrOut = function pointInOrOut (point, index, dir) {
      // assume ray going to + infinity on x plane here just making assumption that it extends 1000 units past whatever the minimum x value is in the boundary
      var p = {
        x0: point[0], y0: point[1], x1: dir, y1: point[1]
      };
      this.ray = p;
      // lets use non-zero winding number rule
      var windingNum = 0;
      var last = { x: Infinity, y: Infinity };

      for (var i = 0; i < index.length; i++) {
        var l = {
          x0: this.coords[index[i]],
          y0: this.coords[index[i] + 1],
          x1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1]],
          y1: this.coords[index[(i + 1) > index.length - 1 ? 0 : i + 1] + 1]
        };
        var inters = intersect(p, l, true);
        if (isFinite(inters.x)) {
          var testCond = Math.round(inters.x * 1000000) === last.x && Math.round(inters.y * 1000000) === last.y;
          if (l.y1 - l.y0 > 0 && !testCond) {
            windingNum++;
          } else if (l.y1 - l.y0 < 0 && !testCond) {
            windingNum--;
          }
          last = { x: Math.round(inters.x * 1000000), y: Math.round(inters.y * 1000000) };
        }
      }
      return Math.abs(windingNum) !== 0
    };

    Boundary.prototype.printPoints = function printPoints (xIndex) {
      var p = [];
      for (var i$1 = 0, list = xIndex; i$1 < list.length; i$1 += 1) {
        var i = list[i$1];

          p.push(this.coords[i], this.coords[i + 1]);
      }
      console.log(p);
    };

    prototypeAccessors.hullCoords.get = function () {
      return this.subset(this.hull)
    };

    Boundary.prototype.subset = function subset (indices) {
      var rv = [];
      for (var i$1 = 0, list = indices; i$1 < list.length; i$1 += 1) {
        var i = list[i$1];

          rv.push(this.coords[i], this.coords[i + 1]);
      }
      return rv
    };

    prototypeAccessors.coords2D.get = function () {
      var newArr = [];
      var arr = this.sortedCoords;
      while (arr.length) { newArr.push(arr.splice(0, 2)); }
      return newArr
    };

    prototypeAccessors.sortedCoords.get = function () {
      var newArr = [];
      for (var i$1 = 0, list = this.hull; i$1 < list.length; i$1 += 1) {
        var i = list[i$1];

          newArr.push(this.coords[i], this.coords[i + 1]);
      }
      return newArr
    };

    Object.defineProperties( Boundary.prototype, prototypeAccessors );

    var BoundaryExtra = /*@__PURE__*/(function (Boundary) {
      function BoundaryExtra (arr, k) {
        if ( k === void 0 ) k = 3;

        Boundary.call(this, arr, k);
        this.cPoints = [];
        this.origCoordsLen = arr.length;
        this.intersectingLineSegs = [];
      }

      if ( Boundary ) BoundaryExtra.__proto__ = Boundary;
      BoundaryExtra.prototype = Object.create( Boundary && Boundary.prototype );
      BoundaryExtra.prototype.constructor = BoundaryExtra;

      var prototypeAccessors = { k: { configurable: true } };

      /**
       * addPoints
       * use final k value from concave boundary for point search in order to
       * @param {Array} parentArr Array of coordinate cloud used to interpolate boundary points
       * @param {Object} constructed delaunator object, need this for edges and triangles
       * @param {Integer} dist Max distance to point to trigger interpolation, only one of two points in line segment has to meet this criteria
       */
      BoundaryExtra.prototype.addPoints = function addPoints (parentArr, delaunator, dist) {
        var this$1 = this;

        this.k = 3;
        var edges = getEdges(delaunator);
        // get all intersecting lines to the hull line seg
        var loop = function ( p ) {
          var h = this$1.subset([this$1.hull[p], this$1.hull[p + 1]]);
          var seg = { x0: h[0], y0: h[1], x1: h[2], y1: h[3] };
          var temp = this$1.getIntersectingLines(seg, edges, parentArr, dist).reverse();
          var ind = sortHeap(temp.map(function (m) { return [m[m.length - 1].x, m[m.length - 1].y]; }).flat(), [].concat( Array(temp.length).keys() ).map(function (i) { return i * 2; }), 'euclid', [seg.x0, seg.y0]);
          this$1.intersectingLineSegs.push(this$1.hull[p]);
          var c = this$1.coords.length;
          this$1.coords = this$1.coords.concat(temp.map(function (m) { return [m[m.length - 1].x, m[m.length - 1].y]; }).flat());
          this$1.intersectingLineSegs = this$1.intersectingLineSegs.concat(ind.map(function (m) { return m + c; }));
        };

        for (var p = 0; p < this$1.hull.length - 1; p++) loop( p );
        this.intersectingLineSegs.push(this.intersectingLineSegs[0]);

        this.hull = this.intersectingLineSegs;
        this.clean(this.hull);
      };

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
      BoundaryExtra.prototype.getIntersectingLines = function getIntersectingLines (lineSeg, indexArr, coords, dist, opt) {
        if ( opt === void 0 ) opt = null;

        var pntAndItsArr = [];
        // iterate over point pairs
        for (var i = 0; i < indexArr.length; i += 2) {
          // first check if distance from either point to perpendicular of line seg is less than dist
          var point = { x: coords[indexArr[i]], y: coords[indexArr[i] + 1] };
          var d = distLineAndPoint(lineSeg, point);
          if (d <= dist) {
            if (i + 1 > indexArr.length - 1) {
              i = -1;
            }
            var coordSeg = { x0: point.x, y0: point.y, x1: coords[indexArr[i + 1]], y1: coords[indexArr[i + 1] + 1] };
            var its = intersect(lineSeg, coordSeg, false);
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
      };

      prototypeAccessors.k.get = function () {
        return Boundary.prototype.k
      };

      prototypeAccessors.k.set = function (v) {
        Boundary.prototype.k = v;
      };

      Object.defineProperties( BoundaryExtra.prototype, prototypeAccessors );

      return BoundaryExtra;
    }(Boundary));

    /**
     * ConstrainoDelaunato
     *
     * @class
     * @classdesc ConstrainoDelaunato
     */
    var ConstrainoDelaunato = function ConstrainoDelaunato (coords, k) {
      var boundaries = [], len = arguments.length - 2;
      while ( len-- > 0 ) boundaries[ len ] = arguments[ len + 2 ];

      // k is the k-nearest neighbor selection
      // if coords are 2D
      if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
        coords = coords.flat();
      } else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
        return
      }

      this.delaunator = new Delaunator(coords);
      this.boundaries = [];
      this.boundedDelaunators = [];

      for (var i = 0, list = boundaries; i < list.length; i += 1) {
        var boundary = list[i];

        if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
          boundary = boundary.flat();
        }
        if (boundary) {
          this.boundaries.push(new BoundaryExtra(boundary, k));
        } else {
          this.boundaries.push(new BoundaryExtra(coords, k));
        }
        this.boundaries[this.boundaries.length - 1].addPoints(coords, this.delaunator, 10);
        this.boundedDelaunators.push(this.setTrianglesInsideBound(this.boundaries[this.boundaries.length - 1]));
      }

      this.boundary = this.boundaries[this.boundaries.length - 1];
      this.boundedDelaunator = this.boundedDelaunators[this.boundedDelaunators.length - 1];
    };

    var prototypeAccessors$1 = { coords2D: { configurable: true },coords: { configurable: true },triangles: { configurable: true },hull: { configurable: true } };

    /**
     * setTrianglesInsideBound
     *
     * Function used to clip coords to inside of boundary or hole
     *
     * @param {BoundaryExtra} boundary boundary extra object
     */
    ConstrainoDelaunato.prototype.setTrianglesInsideBound = function setTrianglesInsideBound (boundary) {
      var coords = [];
      var index = [].concat( this.delaunator.coords.keys() ).filter(function (i) { return i % 2 === 0; });
      var maxX = maximumPointX(this.delaunator.coords, index);
      for (var i$1 = 0, list = index; i$1 < list.length; i$1 += 1) {
        var e = list[i$1];

          var point = { x: this.delaunator.coords[e], y: this.delaunator.coords[e + 1] };
        // if (point.x === 59 && point.y === 80) {
        // console.log(boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10))
        // }
        if (boundary.pointInOrOut([point.x, point.y], boundary.hull, maxX.x + 10)) {
          coords.push(point.x, point.y);
        }
      }

      coords = coords.concat(boundary.subset(boundary.hull));
      var rv = new Delaunator(coords);
      var t = [];
      for (var e$1 = 0; e$1 < rv.triangles.length / 3; e$1++) {
        var edgeIndex = e$1 * 3;

        var xCoord = 0;
        var yCoord = 0;
        for (var r = 0; r < 3; r++) {
          xCoord += rv.coords[2 * rv.triangles[r + edgeIndex]];
          yCoord += rv.coords[2 * rv.triangles[r + edgeIndex] + 1];
        }
        var point$1 = { x: xCoord / 3, y: yCoord / 3 };
        if (boundary.pointInOrOut([point$1.x, point$1.y], boundary.hull, maxX.x + 10)) {
          t.push(rv.triangles[edgeIndex], rv.triangles[edgeIndex + 1], rv.triangles[edgeIndex + 2]);
        }
      }
      rv.triangles = new rv.triangles.constructor(t);

      return rv
    };

    /**
     * update
     *
     * @param {Array} point x and y coord of point to add the delaunator object
     */
    ConstrainoDelaunato.prototype.update = function update (point) {
      var c = this.coords;
      for (var i = 0, list = point.flat(); i < list.length; i += 1) {
        var p = list[i];

          c.push(p);
      }
      this.delaunator = new Delaunator(c);
    };

    /**
     * coords2D
     *
     * @returns {Array} 2D coordinate array
     */
    prototypeAccessors$1.coords2D.get = function () {
      var c2D = [];
      var c1D = this.coords;
      for (var i = 0; i < c1D.length; i += 2) {
        c2D.push([c1D[i], c1D[i + 1]]);
      }
      return c2D
    };

    /**
     * coords
     *
     * @returns {Array} 1D coordinate array
     */
    prototypeAccessors$1.coords.get = function () {
      return this.delaunator.coords
    };

    /**
     * triangles
     *
     * @returns {Array} Index array of delaunator triangles
     */
    prototypeAccessors$1.triangles.get = function () {
      return this.delaunator.triangles
    };

    /**
     * hull
     *
     * @returns {Array} Array of hull indices
     */
    prototypeAccessors$1.hull.get = function () {
      return this.delaunator.hull
    };

    Object.defineProperties( ConstrainoDelaunato.prototype, prototypeAccessors$1 );

    return ConstrainoDelaunato;

})));
