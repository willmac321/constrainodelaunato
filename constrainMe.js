import Delaunator from 'delaunator';

const test = [10,6,3,4,7,1,2,5];

export default class ConstrainoDelaunato{
	constructor(coords, boundary, k) {
		//k is the k-nearest neighbor selection
		// if coords are 2D
		if (coords && Array.isArray(coords[0]) && coords[0].length === 2) {
			coords = coords.flat();
		} else if (coords && Array.isArray(coords[0]) && coords[0].length !== 2) {
			return
		}
		if (boundary && Array.isArray(boundary[0]) && boundary[0].length === 2) {
			boundary = boundary.flat();
		}
		if(boundary) {
	//		boundary = this.sortHeap(boundary);
			this.sortHeap(test, 1)
			this.boundary = this.makaThaEnvelope(boundary);
			coords = coords.concat(boundary);
		}
		this.bound = boundary;
		this.boundary = new Delaunator(boundary);
		this.delaunator = new Delaunator(coords);
//		this.pointInOrOut([1,1]);
	}

	pointInOrOut(point) {
		//lets use non-zero winding number rule
		let windingNum = 0
		for (let h of this.boundary.hull) {
			console.log(h);
		}
//		console.log(this.coords);
//		for (let e = 0; e < this.delaunator.triangles.length; e++) {
//			if (e > this.delaunator.halfedges[e]) {
//				let p = points[this.delaunator.triangles[e]];
//				let q = points[this.delaunator.triangles[nextHalfedge(e)]];
//			}
//		}
	}

	makaThaEnvelope(arr, k) {
//k nearest neighbor babbbbyyyy
//https://towardsdatascience.com/the-concave-hull-c649795c0f0f
		

	}

	sortHeap(arr, dim) {
		let newArr = [];
		let index = 0;
		let minY = Infinity;
		let minX = Infinity;
		//convert point arr to 2d -> easier for me to get my head around sorting
		if (dim > 1) {
			while(arr.length) newArr.push(arr.splice(0, dim));
		}
		else{
			newArr = arr.slice();
		}

		for (let p = 0; p < newArr.length; p++) {
			//TODO remove
			if (newArr[p] < minY) {
				minY = newArr[p];
				index = p;
			}

//			if (newArr[p][1] < minY) {
//				minX = newArr[p][0];
//				minY = newArr[p][1];
//				index = p;
//			}
//			else if (newArr[p][1] <= minY && newArr[p][0] <= minX) {
//				minX = newArr[p][0];
//				minY = newArr[p][1];
//				index = p;
//			}
		}
		console.log(minY, newArr.slice());
//		builtInSort([minX, minY], newArr);
		//		TODO change back
//		heapSort([minX, minY], newArr, newArr.length);
		heapSort(minX, newArr, newArr.length);

		return newArr.flat();
	}

	update(point) {
		let c = this.coords;
		for (let p of point.flat()) {
			c.push(p);
		}
		this.delaunator = new Delaunator(c);
	}

	get coords2D() {
		let c2D = []
		let c1D = this.coords;
		for (let i = 0; i < c1D.length; i += 2) {
			c2D.push(this.c1D[i], this.c1D[i + 1]);
		}
		return c2D;
	}

	get coords() {
		return this.delaunator.coords;
	}

	get triangles() {
		return this.delaunator.triangles;
	}

	get hull() {
		return this.delaunator.hull;
	}
}

function nextHalfEdge(e) {
	return (e % 3 === 2) ? e - 2 : e + 1;
}

function dotProduct(a, b) {
	let p = {x: a[0], y: a[1]};
	let o = {x: b[0], y: b[1]};
	return p.x*o.x + p.y*o.y;
}

function manhattenDist(a, b) {
	let p = {x: a[0], y: a[1]};
	let o = {x: b[0], y: b[1]};
	return Math.abs(p.x - o.x) + Math.abs(p.y - o.y);
}

function builtInSort(point, arr) {
	return arr.sort((a, b) => {
//		return manhattenDist(point, a) - manhattenDist(point, b)
		return (dotProduct(point, a) - dotProduct(point, b));
	});
}

//heap sort 2d array by angle
function heapSort(point, a, count, p) {
	let func;
	let prom = new Promise((res, rej) => {
		if (p === 'dist'){
			func = manhattenDist;
		} else if (!Array.isArray(a[0])) {
			console.log('here');
			func = (a, b) => a - b;
		} else {
			func = dotProduct;
		}
		return res();
	});
	
	prom.then((res)=> {
		heapify(point, a, count, func);
		return console.log(point, a, count);
	}).then((r)=> {
		let end = count - 1;
		while (end > 0) {
			swap(a, end, 0);
			end--;
			siftDown(point, a, 0, end, func);
		}
	}).then(() => {
		for (let r of a) {
			console.log(r, func(r, point));
		}
	}).catch((e) => console.log('uhoh'));
}

function heapify(point, a, count, func) {
	let par = (i) => Math.floor((i - 1) / 2);
	let start = par(count - 1);
	while (start >= 0) {
		siftDown(point, a, start, count - 1, func);
		start--;
	}
}

function siftDown(point, a, start, end, func) {
	let root = start;
	let left = (i) => 2 * i + 1;
	let right = (i) => 2 * i + 2;

	while (left(root) <= end) {
		let child = left(root);
		let s = root;
		let dot = (i) => { 
//			return dotProduct(a[i], point);
			let t = func(a[i], point);
//			console.log(a[i], point, t);
			return t;
		}
		let dotOld = dot(s); 

		if (dotOld < dot(child)) {
			s = child;
		}
		if (child + 1 <= end && dotOld < dot(child + 1)) {
			s = child + 1;
		}
		if (s === root) {
			return;
		} else {
//			console.log(`swap ${a[root]} ${a[s]}`);
			swap(a, root, s);
			root = s;
		}
	}
}

function swap(a, i, j) {
	let t = a[i];
	a[i] = a[j];
	a[j] = t;
}

