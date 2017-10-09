
template<class T>
inline std::string toString(const T& t){
	std::ostringstream os;
	os << t;
  	return os.str();
}

template<class T>
inline T fromString(const std::string& s){
	T t;
	std::istringstream is(s);
	is >> t;
	return t;
}

std::vector<double> DeltaR(std::vector<double> a, std::vector<double> b, std::vector<double> l){ // distance
	std::vector<double> r (3,0);
        for (int z=0; z<3; z++)	r.at(z) = a.at(z) - b.at(z); // boundary correction
        if (r[0] > l[0]*0.5)  r[0] -= l[0];
        if (r[1] > l[1]*0.5)  r[1] -= l[1];
        if (r[2] > l[2]*0.5)  r[2] -= l[2];

        if (r[0] < -l[0]*0.5)  r[0] += l[0];
        if (r[1] < -l[1]*0.5)  r[1] += l[1];
        if (r[2] < -l[2]*0.5)  r[2] += l[2];
	return r;
}

double scal_p (std::vector<double> a1, std::vector<double> a2){ // scalar product
	double r=0;
	for(int v=0; v<3; v++) r += a1[v]*a2[v];
	return r;
}


std::vector<double> kreuz(std::vector<double> a, std::vector<double> b, std::vector<double> c){ // distance
    std::vector<double> oo (3,0);
    oo[0] = a[1]*b[2] - a[2]*b[1];
    oo[1] = a[2]*b[0] - a[0]*b[2];
    oo[2] = a[0]*b[1] - a[1]*b[0];

    double o = scal_p(oo,c);
    if (o > 0){
        oo[0] *=-1;
        oo[1] *=-1;
        oo[2] *=-1;
    }
    return oo;
}
