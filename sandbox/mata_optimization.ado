version 12

mata: mata clear

mata:



	void eval0(todo, p, v, g, H){
		v = -(p[1]-4)^2 - (p[2]-5)^2 - (p[2]-p[1])^4;
	}

	S = optimize_init()
	optimize_init_evaluator(S, &eval0())
	optimize_init_params(S, (0,0))
	p = optimize(S)
	p

	void eval1(todo, p, v, g, H){
		v = -(p[1]-4)^2 - (p[2]-5)^2 - (p[2]-p[1])^4;
		if(todo>=1){
			g[1] = -2 *(p[1]-4) - 4* (p[2]-p[1])^3 *(-1);
			g[2] = -2 *(p[2]-5) - 4* (p[2]-p[1])^3 *(1);
		}
	}

	S = optimize_init()
	optimize_init_evaluator(S, &eval1())
	optimize_init_evaluatortype(S, "d1")
	optimize_init_params(S, (0,0))
	p = optimize(S)
	p



end
