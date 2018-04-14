version 13

drop _all
mata: mata clear

cd "C:\Users\David\YandexDisk\hsework\gauss-mata\sandbox"

mata
// a couple of programs that allow to write .csv (without headers)
function writeMatrix(filename, data) {
	unlink(filename)
	fh = fopen(filename, "w")
	for (i = 1; i <= rows(data); i++) {
		fput(fh, invtokens(strofreal(data[i,]), ","))
	}
	fclose(fh)
}

function readMatrix(filename) {
	fh = fopen(filename, "r")
	tok = tokeninit(",")
	if ((line=fget(fh))!=J(0,0,"")) {
		tokenset(tok, line)
		result = strtoreal(tokengetall(tok))
	}
	while ((line=fget(fh))!=J(0,0,"")) {
		tokenset(tok, line)
		result = result \ strtoreal(tokengetall(tok))
	}
	fclose(fh)
	return(result)
}

x = rnormal(30, 1, 5,2)
y = 15 :+ 0.2*x + rnormal(30,1,0,2)

writeMatrix("test.csv", (x,y))

readMatrix("test.csv")

end
