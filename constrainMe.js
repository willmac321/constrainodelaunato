import Delaunator from 'delaunator';
import {heapSort} from './helpers';

const test = [10,6,3,4,7,1,2,5];

class Boundary{
	constructor(arr) {
		this.coords = arr.slice();
		console.log(this.coords);
		this.indices = [];
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
			this.boundary = new Boundary(sortHeap(boundary, 2));
			sortHeap(test, 1)
			coords = coords.concat(this.boundary.coords);
		}
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

function sortHeap(arr, dim) {
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
		console.log(minX, newArr.slice());
//		builtInSort([minX, minY], newArr);
		heapSort([minX, minY], newArr, newArr.length);
		console.log(minX, newArr.slice());

		return newArr.flat();
}

