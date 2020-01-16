import Delaunator from 'delaunator';
import {heapSort} from './helpers';

const test = [10,6,3,4,7,1,2,5];

class Boundary{
	constructor(arr) {
		this.coords = arr.slice();
		this.indices = [];
		this.center = this.calcCenter();
		this.coords = sortHeap(this.coords, 2, 'polar', [this.center.x, this.center.y]);
		this.minY = minimumPointY(this.coords2D);
		this.minX = minimumPointX(this.coords2D);
		this.ray;
		//center of test from html is not inside boundary
		// this point is though
		
		this.pointInOrOut([this.center.x, this.center.y]);
		this.pointInOrOut([this.minX.x + 1000, this.minX.y]);
		console.log(this.pointInOrOut([180, 100]));
	}

	//concave center point
	calcCenter() {
		let p = {x: 0, y: 0};

		for (let i = 0; i < this.coords.length; i += 2) {
			p.x += this.coords[i];
			p.y += this.coords[i + 1];
		}
		p.x = p.x / (this.coords.length / 2)
		p.y = p.y / (this.coords.length / 2)
		return p;
	}

	pointInOrOut(point) {

		//assume ray going to + infinity on x plane
		let p = {x0: point[0], y0: point[1], x1: this.minX.x - 1000, y1: point[1]};
		this.ray = p;
		console.log(this.ray);
		//lets use non-zero winding number rule
		let windingNum = 0

		for (let i = 0; i < this.coords.length; i += 2) {
			let l = {
				x0: this.coords[i], 
				y0: this.coords[i + 1], 
				x1: this.coords[(i + 2 ) > this.coords.length - 1 ? 0 : i + 2],
				y1: this.coords[(i + 3 ) > this.coords.length - 1 ? 1 : i + 3]
			}
//			console.log(l, p);
			let intersect = this.intersect(p, l);
			if (isFinite(intersect.x)) {
				if(l.y1 - l.y0 > 0) {
					windingNum++;
				} else if(l.y1 - l.y0 < 0) {
					windingNum--;
				}
			} 
		}
		console.log(windingNum);
		return windingNum !== 0;
	}

	intersect(p, l) {
		let den = ((l.y1 - l.y0) * (p.x1 - p.x0)) - ((l.x1 - l.x0) * (p.y1 - p.y0));
		if (den === 0) {
			return {x:Infinity, y:Infinity};
		}

		let a = p.y0 - l.y0;
		let b = p.x0 - l.x0;

		let num1 = ((l.x1 - l.x0) * a) - ((l.y1 - l.y0) * b);
		let num2 = ((p.x1 - p.x0) * a) - ((p.y1 - p.y0) * b);

		a = num1 / den;
		b = num2 / den;

		let rv = {
			x: p.x0 + (a * (p.x1 - p.x0)),
			y: p.y0 + (a * (p.y1 - p.y0))
		};

		let t = {a:false, b:false};

//		if (p.y1 === rv.y) {
//			console.log(a, b);
//		}
//
		if (a >= 0 && a < 1) {
			t.a = true;	
		}
		if (b >= 0 && b < 1) {
			t.b = true;	
		} 
		if(t.a && t.b) {
			return rv;
		} else {
			return {x:Infinity, y:Infinity};
		}
	}


	makaThaEnvelope(arr, k) {
//k nearest neighbor babbbbyyyy
//https://towardsdatascience.com/the-concave-hull-c649795c0f0f

	}

	get coords2D() {
		let newArr = [];
		let arr = this.coords.slice();
		while(arr.length) newArr.push(arr.splice(0, 2));
		return newArr;
	}
}

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
			this.boundary = new Boundary(boundary);
//			sortHeap(test, 1)
			coords = coords.concat(this.boundary.coords);
		}
		this.delaunator = new Delaunator(coords);
//		this.pointInOrOut([1,1]);
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
			c2D.push([c1D[i], c1D[i + 1]]);
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

	get bound() {
		return this.boundary.coords;
	}
}

function sortHeap(arr, dim, criteria, centerPoint) {
		let newArr = [];
		//convert point arr to 2d -> easier for me to get my head around sorting
		if (dim > 1) {
			while(arr.length) newArr.push(arr.splice(0, dim));
		}
		else{
			newArr = arr.slice();
		}
		let minPoint = minimumPointY(newArr);
		console.log(minPoint.x, minPoint.y, newArr.slice());
//		builtInSort([minX, minY], newArr);
		heapSort([minPoint.x, minPoint.y], newArr, newArr.length, criteria, centerPoint);
		console.log(minPoint.x, minPoint.y, newArr.slice());

		return newArr.flat();
}

function minimumPointY(newArr) {
		let index = 0;
		let minY = Infinity;
		let minX = Infinity;
		if (Array.isArray(newArr[0])) {
			for (let p = 0; p < newArr.length; p++) {
				if (newArr[p][1] < minY) {
					minX = newArr[p][0];
					minY = newArr[p][1];
					index = p;
				}
				else if (newArr[p][1] <= minY && newArr[p][0] <= minX) {
					minX = newArr[p][0];
					minY = newArr[p][1];
					index = p;
				}
			}
		} else {
			for (let p = 0; p < newArr.length; p++) {
				if (newArr[p] < minX) {
					minX = newArr[p];
					index = p;
				}
			}
		}
	return {x: minX, y:minY};
}

function minimumPointX(newArr) {
		let index = 0;
		let minY = Infinity;
		let minX = Infinity;
		if (Array.isArray(newArr[0])) {
			for (let p = 0; p < newArr.length; p++) {
				if (newArr[p][0] < minX) {
					minX = newArr[p][0];
					minY = newArr[p][1];
					index = p;
				}
				else if (newArr[p][1] <= minY && newArr[p][0] <= minX) {
					minX = newArr[p][0];
					minY = newArr[p][1];
					index = p;
				}
			}
		} else {
			for (let p = 0; p < newArr.length; p++) {
				if (newArr[p] < minX) {
					minX = newArr[p];
					index = p;
				}
			}
		}
	return {x: minX, y: minY};
}
