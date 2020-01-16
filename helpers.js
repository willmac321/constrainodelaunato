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
export function heapSort(point, a, count, p) {
	let func;
	if (p === 'dist'){
		func = manhattenDist;
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

