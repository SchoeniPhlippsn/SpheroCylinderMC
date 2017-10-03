void Init(){


    std::ifstream status("Save/status.txt");

    Config.rhoN = rho0/Config.Vsys;
    Config.Vbox = N/Config.rhoN;
        
    Config.l[0] =pow(Config.Vbox,1.0/3.0);
    Config.l[1] = Config.l[0];
    Config.l[2] = Config.l[0];

	Config.W.resize(3,0); 
	Config.w.resize(3,0); 

    Config.W[0] = (int)(Config.l[0]/rlist); 
    Config.W[1] = (int)(Config.l[1]/rlist); 
    Config.W[2] = (int)(Config.l[2]/rlist); 

    if(status) status >> restart;
    else{
        std::cerr << "There is no Save/status.txt file! I don't know if I have to restart." << std::endl;
        exit(-1);
    }

    if(restart){ 
        Config.read("Save/Config.dat", Rdist, rdist);
        initCompress=false;
    }else{
        initCompress=false;

        Config.rhoV = 0.01;
        Config.rhoN = Config.rhoV/Config.Vsys;
        Config.Vbox = N/Config.rhoN;
            
        Config.l[0] =pow(Config.Vbox,1.0/3.0);
        Config.l[1] = Config.l[0];
        Config.l[2] = Config.l[0];
        
        std::cout << "System size\tlx=" << Config.l[0] << "\tly=" << Config.l[1] << "\tlz=" << Config.l[2] << std::endl;



        int grid = pow(N,1.0/3.0);
        grid++; 
        
        int i_min = 1;
        int i = 1;
        int j = 1;
        int k = 1;
        int bk = 1;
        double stepi = Config.l[0]/grid;
        double stepj = Config.l[1]/grid;
        double stepk = Config.l[2]/grid;

        for (int ii=0; ii < N; ii++){			
            bool inside =false;
            double rmin = 1000;
            int h=0;
            class teilchen newpart;
            newpart.pos.push_back( i*stepi);
            newpart.pos.push_back( j*stepj);
            newpart.pos.push_back( k*stepk);
            if( k % 2 == 0 ){ 
                newpart.pos[0] += 0.5*stepi; 
                newpart.pos[1] += 0.5*stepj; 
            }
            if(ii < Config.Nc ){
                newpart.rdist = rdist; 
                newpart.Rdist =-Rdist; 
                newpart.dist = rdist-Rdist;
                newpart.R = Rsp;
                newpart.r = rsp;

                newpart.ori.push_back(0); 
                newpart.ori.push_back(0); 
                newpart.ori.push_back(bk); 
            
                for( int v=0; v< Config.part.size(); v++){
                    if(v < Config.Nc ){
                        if(overlapSP(newpart,Config.part[v],Config.l)){
                            inside = true;
                            break;   
                        }
                    }else{
                        if(overlapSPSPH(newpart,Config.part[v],Config.l)){
                            inside = true;
                            break;   
                        }
                    }
                }	
            }else{
                newpart.r = rsphere;
                for( int v=0; v< Config.part.size(); v++){
                    if(v < Config.Nc ){
                        if(overlapSPSPH(Config.part[v],newpart,Config.l)){
                            inside = true;
                            break;   
                        }
                    }else{
                        if(overlapSPH(newpart,Config.part[v],Config.l) ){
                            inside = true;
                            break;   
                        }
                    }
                }	
            }
            if(inside){
                std::cerr << "Failed to place particle! Intersection! " << grid << " " << i << "  " << j << " " << k << std::endl;
                exit(1);
            }


            Config.part.push_back(newpart);

            i+=2;
            if (i > grid && ii != N-1){
                j++;
                i=i_min;
                if (j > grid){
                    k++;
                    bk*=-1;
                    j=1;
                    if(k>grid){
                        if( i_min == 2 ){
                            std::cerr << "There is a problem with the placement of the particles!" << std::endl;
                            exit(1);
                        }
                        std::cout << ii << std::endl;
                        j = 1;
                        k = 1;
                        i_min =2;
                        i = 2;
                    }
                }
            }
        }
    }

    RenewList();

}
