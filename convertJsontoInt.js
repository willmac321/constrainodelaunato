const fs = require('fs');

function dothething(jso) {

	let raw  = fs.readFileSync(jso);
	let obj = JSON.parse(raw);

	for (let i = 0, j = 0; i < obj.length; i ++) {
		obj[i][j] = Math.floor(obj[i][j]);
		j++;	
		obj[i][j] = Math.floor(obj[i][j]);
		j--;	
		if (i >= 1 && obj[i][0] === obj[i - 1][0] && obj[i][1] === obj[i - 1][1]) {
			obj.splice(i, 1);
			i--;
		}
	}
	
	fs.writeFileSync('temp' + process.argv[2], JSON.stringify(obj));

}

dothething(process.argv[2]);
