var c; 

function nextHalfEdge(e) {
	return (e % 3 === 2) ? e - 2 : e + 1;
}

function dotProduct(a, b) {
	let p = {x: a[0], y: a[1]};
	let o = {x: b[0], y: b[1]};
	return p.x*o.x + p.y*o.y;
}

function dotPolar(a, b) {
	// b is basis point
	// put those bad boys in order ccw around some centroid point that globally declared
	let o = {x: c.x - a[0], y: c.y - a[1]};
	let p = {x: c.x - b[0], y: c.y - b[1]};
	let magO = Math.sqrt(Math.pow(o.x, 2) + Math.pow(o.y, 2));
	let magP = Math.sqrt(Math.pow(p.x, 2) + Math.pow(p.y, 2));
	let theta = Math.acos((p.x*o.x + p.y*o.y) / (magO * magP)) * 180 / Math.PI;
	let det = (o.x * p.y) - (p.x * o.y);
	theta = isNaN(theta) ? 0 : theta;
	theta = det > 0 ? 360 - theta : theta;
	return theta;
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
export function heapSort(point, a, count, p, center) {
	let func;
	if (p === 'dist'){
		func = manhattenDist;
	} else if (p === 'polar') {
		func = dotPolar;
		c = {x: center[0], y: center[1]};
	} else if (!Array.isArray(a[0])) {
		point = point[0];
		func = (a, b) => a - b;
	} else {
		func = dotProduct;
	}

	heapify(point, a, count, func);

	let end = count - 1;
	while (end > 0) {
		swap(a, end, 0);
		end--;
		siftDown(point, a, 0, end, func);
	}
//	for (let r of a) {
//		console.log(r, func(r, point));
//	}
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

	while (left(root) <= end) {
		let child = left(root);
		let s = root;
		let dot = (i) => { 
			let t = func(a[i], point);
			return t;
		};

		if (dot(s) < dot(child)) {
			s = child;
		}
		if (child + 1 <= end && dot(s) < dot(child + 1)) {
			s = child + 1;
		}
		if (s === root) {
			return;
		} else {
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

