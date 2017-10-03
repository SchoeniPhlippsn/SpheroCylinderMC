#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
	
boost::random::mt19937 gen;	// random number generator [0,1]
boost::random::uniform_01<> uni;


#include "Headers/Headers.h"



int main(int argc, char** argv){

	Config.Nc=2; //# of spherocylinder
	Config.Ns=498; //# of spheres
	seed=42;
	rho0=0.45; //density
	pos_lambda = 0.02;
	ori_lambda = 0.02;
	vproc = 0.99;
    Config.l.push_back(0);
    Config.l.push_back(0);
    Config.l.push_back(0);
    Rdist = -0.884;
    rdist = 1.233;
    Rsp = 0.616;
    rsp = 0.267;
    Vratio = 0.1;
    finstep = 500000;
    savestep = 500;
    backupstep = 100;

	for (int argi=1; argi<argc; argi++){
        std::string argvi = argv[argi];
        int next = argc-argi-1;
        if (argvi=="-Ns" && next>0) Config.Ns      = fromString<long>(argv[++argi]);	
        else if (argvi=="-Nc" && next>0) Config.Nc      = fromString<long>(argv[++argi]);	
		else if (argvi=="-s" && next>0) seed  = fromString<long>(argv[++argi]);
		else if (argvi=="-fs" && next>0) finstep  = fromString<long>(argv[++argi]);
		else if (argvi=="-sa" && next>0) savestep  = fromString<long>(argv[++argi]);
		else if (argvi=="-v" && next>0) vproc = fromString<double>(argv[++argi]);
		else if (argvi=="-Vr" && next>0) Vratio = fromString<double>(argv[++argi]);
		else if (argvi=="-Rd" && next>0) Rdist = fromString<double>(argv[++argi]);
		else if (argvi=="-rd" && next>0) rdist = fromString<double>(argv[++argi]);
		else if (argvi=="-Rsp" && next>0) Rsp = fromString<double>(argv[++argi]);
		else if (argvi=="-rsp" && next>0) rsp = fromString<double>(argv[++argi]);
		else if (argvi=="-r" && next>0) rho0 = fromString<double>(argv[++argi]);
		else if (argvi=="-pl" && next>0) pos_lambda = fromString<double>(argv[++argi]);
		else if (argvi=="-ol" && next>0) ori_lambda = fromString<double>(argv[++argi]);
		else throw std::runtime_error("no such argument: "+argvi);
  	}

    RHO = toString(rho0);
    NC = toString(Config.Nc);
    NS = toString(Config.Ns);
    VR = toString(Vratio);

    height = rdist-Rdist;
    Rsp_cube= Rsp*Rsp*Rsp;
    rsp_cube= rsp*rsp*rsp;

    Config.Vsys = (rsp_cube+Rsp_cube)*2.0/3.0*M_PI;

    Config.Vsys += M_PI*height/(3.0*(Rsp-rsp))*(Rsp_cube-rsp_cube);

    std::cout << "V_c="<< Config.Vsys << std::endl;

    rsphere = pow(3.0/4.0*Vratio*Config.Vsys/M_PI,1.0/3.0);

    std::cout << rsphere << std::endl;

    Config.Vsys = (Config.Nc + Config.Ns*Vratio)*Config.Vsys;

    N = Config.Nc + Config.Ns;

    Config.Vsys /= N;

    rlist = 2*rsphere + 0.1;
    rlist_2 = rlist*rlist;  
    maxpos = 0.2*rlist;

	gen.seed(seed);

    Config.step = 0;
    std::cout << "Starting Initialisation of the system!" << std::endl; 
	Init();
    std::cout << std::endl;
    std::cout << "Particles are placed successfully!" << std::endl;

    Compression();

    Config.step++;

    acc = 0;   
    while(Config.step <= finstep){

        Move_step();
        if(Config.step % savestep == 0){          
            std::string file_name = "Results/Config_Nc" + NC + "_Ns"+ NS + "_Vr" + VR + "_rho" + RHO + "_step"+ toString(Config.step) + ".dat";
            Config.write(file_name,0);
        }
        if(Config.step % backupstep == 0){          
            Config.write("Save/Config.dat",1);
        }

        if(Config.step % 100 == 0){          
            acceptance = static_cast<double>(acc)/(100*N);
            std::cout << Config.step << " " << acceptance << std::endl;
            if (acceptance > 0.55){
                if( pos_lambda < maxpos) pos_lambda += 0.1;
                else pos_lambda = maxpos;
                if( ori_lambda < maxpos) ori_lambda += 0.02;
                else ori_lambda = maxpos;
            }
            if (acceptance < 0.45){
                if( pos_lambda < 0.01 ) pos_lambda *= 0.5;
                else pos_lambda -= 0.01;
                if( ori_lambda > 0.01 ) ori_lambda -= 0.01;
            }
            acc=0;
        }
        Config.step++;
    }
    Config.write("Results/finalConfig.dat",1);

    std::ofstream File ("Save/status.txt");

    File << 0 << std::endl;

    File.close();
}
