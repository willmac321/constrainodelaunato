import Delaunator from 'delaunator';

export default class ConstrainoDelaunato{
	constructor(coords, boundary) {
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
			coords = coords.concat(boundary);
		}
		this.boundary = boundary;
		console.log(coords);
		this.delaunator = new Delaunator(coords);
	}

	pointInOrOut(point) {
		//lets use non-zero winding number rule
		let windingNum = 0


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
