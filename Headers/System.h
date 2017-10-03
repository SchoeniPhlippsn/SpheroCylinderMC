class teilchen{
    public:
        
    std::vector<double> pos;       	// position of the particle
    std::vector<double> ori;      	// orientation of the particle
    double R;
    double r;
    double dist;
    double rdist;
    double Rdist;
    double a;
    double c;
	int cell;

        private:
};


class system{
public:
    
    std::vector<teilchen> part;       
	std::vector<double> l;
	double rhoN;
	double rhoV;
	double Vbox;
	double Vsys;
    int Nc;
    int Ns;
    std::vector<int> head;
    std::vector<int> link;
    std::vector<int> W;
    std::vector<double> w;
    int step;
    void write(std::string, bool);
    void read(std::string,double,double);
	                                                     
	private:
    int v,vv;
};

void system::write (std::string file, bool save){

    std::ofstream oFile (file.c_str() );
    if(!save){
        oFile << 2*Nc + Ns << std::endl;
        oFile << 0 << " " << l[0] << std::endl;
        oFile << 0 << " " << l[1] << std::endl;
        oFile << 0 << " " << l[2] << std::endl;
    }else{
        oFile << Nc << " " << Ns << std::endl;
        oFile << l[0] << " " << l[1] << " " << l[2] << std::endl;
        oFile << step << std::endl;

        std::ofstream ooFile("Save/status.txt");
        ooFile << 1 << std::endl;
        ooFile.close();
    }

    for( v = 0 ; v < Nc; v++){

        if(!save){ 
            oFile << "SPHERE ";
            oFile <<  part[v].pos[0] + part[v].rdist*part[v].ori[0] << " " << part[v].pos[1] + part[v].rdist*part[v].ori[1] << " " << part[v].pos[2] + part[v].rdist*part[v].ori[2] << " ";
            oFile << 2*part[v].r << " 0 0 0 " << 2*part[v].r << " 0 0 0 " << 2*part[v].r << " " << 1 <<std::endl;
            oFile << "SPHERE ";
            oFile <<  part[v].pos[0] - part[v].Rdist*part[v].ori[0] << " " << part[v].pos[1] - part[v].Rdist*part[v].ori[1] << " " << part[v].pos[2] - part[v].Rdist*part[v].ori[2] << " ";
            oFile << 2*part[v].R << " 0 0 0 " << 2*part[v].R << " 0 0 0 " << 2*part[v].R << " " << 2 << std::endl;
        }else{
            oFile <<  part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
            oFile <<  part[v].ori[0] << " " << part[v].ori[1] << " " << part[v].ori[2] << " "; 
            oFile <<  part[v].R << " " << part[v].r << " " << part[v].Rdist << " " << part[v].rdist << std::endl; 
        }
    }
    for( v = Nc ; v < part.size(); v++){

        if(!save){ 
            oFile << "SPHERE ";
            oFile <<  part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " ";
            oFile << 2*part[v].r << " 0 0 0 " << 2*part[v].r << " 0 0 0 " << 2*part[v].r << " " << v << std::endl;
        }else{
            oFile << part[v].pos[0] << " " << part[v].pos[1] << " " << part[v].pos[2] << " "; 
            oFile << part[v].r <<  std::endl; 
        }
    }
    oFile.close();
}

void system::read (std::string file, double Rdist, double rdist ){

    std::ifstream iFile (file.c_str() );

    if(!iFile){
        std::cerr << "Can not read " << file << "!" << std::endl;
        exit(-1);
    }
    iFile >> Nc >> Ns;
    part.resize(Nc+Ns);
    l.resize(3);
    
    iFile >> l[0] >> l[1] >> l[2];

    iFile >> step;

    Vbox=l[0]*l[1]*l[2];

    rhoN = part.size()/Vbox;
    rhoV = Vsys*rhoN;

    for( v = 0 ; v < Nc; v++){
        part[v].pos.resize(3);
        part[v].ori.resize(3);

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2] >> part[v].ori[0] >> part[v].ori[1] >> part[v].ori[2] >> part[v].R >> part[v].r >> part[v].Rdist >> part[v].rdist; 

        part[v].dist=part[v].Rdist+part[v].rdist;

        part[v].c = (part[v].R-part[v].r)/part[v].dist;
        part[v].a = part[v].R - part[v].Rdist*part[v].c;
    }
    for( v = Nc ; v < part.size(); v++){
        part[v].pos.resize(3);
        part[v].ori.resize(3);

        iFile >> part[v].pos[0] >> part[v].pos[1] >> part[v].pos[2] >> part[v].r; 
    }
    iFile.close();
}
