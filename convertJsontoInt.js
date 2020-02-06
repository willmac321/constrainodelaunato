const fs = require('fs')

function dothething (jso) {
  const raw = fs.readFileSync(jso)
  const obj = JSON.parse(raw)

  for (let i = 0, j = 0; i < obj.length; i++) {
    obj[i][j] = Math.floor(obj[i][j])
    j++
    obj[i][j] = Math.floor(obj[i][j])
    j--
    if (i >= 1 && obj[i][0] === obj[i - 1][0] && obj[i][1] === obj[i - 1][1]) {
      obj.splice(i, 1)
      i--
    }
  }

  fs.writeFileSync('temp' + process.argv[2], JSON.stringify(obj))
  console.log(`file saved as temp${jso}`)
}

function dothething2 (jso) {
  const raw = fs.readFileSync(jso)
  const obj = JSON.parse(raw)
  const o = []

  for (let i = 0; i < obj.length; i += 2) {
    o.push([obj[i], obj[i + 1]])
  }

  fs.writeFileSync('temp' + process.argv[2], JSON.stringify(o))
  console.log(`file saved as temp${jso}`)
}

dothething2(process.argv[2])
