void Move_step(){
    for(int v = 0;  v < N; v++){
//        std::cout << "Step: " << v << std::endl; 
      
        randP = uni(gen)*N;
        if (randP==N){
            v--;
            continue;
        }
       
        MovedParticle = Config.part[randP];         

        bool inside = false;

        if( randP < Config.Nc ){

            if( 2*uni(gen) < 1 ){
                MovedParticle.pos[0] = MovedParticle.pos[0] + (uni(gen)-0.5)*pos_lambda;
                MovedParticle.pos[1] = MovedParticle.pos[1] + (uni(gen)-0.5)*pos_lambda;
                MovedParticle.pos[2] = MovedParticle.pos[2] + (uni(gen)-0.5)*pos_lambda;


            }else{
                double phi = (uni(gen)-0.5)*ori_lambda;
                double cosphi = cos(phi);
                double sinphi = sin(phi);

                int Case = 3*uni(gen);
                switch(Case){
                case 0:
                    MovedParticle.ori[0] = cosphi*MovedParticle.ori[0] + sinphi*MovedParticle.ori[1];
                    MovedParticle.ori[1] = -sinphi*MovedParticle.ori[0] + cosphi*MovedParticle.ori[1];
                    break;
                case 1:
                    MovedParticle.ori[2] = cosphi*MovedParticle.ori[2] + sinphi*MovedParticle.ori[0];
                    MovedParticle.ori[0] = -sinphi*MovedParticle.ori[2] + cosphi*MovedParticle.ori[0];
                    break;
                case 2:
                    MovedParticle.ori[1] = cosphi*MovedParticle.ori[1] + sinphi*MovedParticle.ori[2];
                    MovedParticle.ori[2] = -sinphi*MovedParticle.ori[1] + cosphi*MovedParticle.ori[2];
                    break;
                }

                double norm = 1/sqrt(MovedParticle.ori[0]*MovedParticle.ori[0] + MovedParticle.ori[1]*MovedParticle.ori[1] + MovedParticle.ori[2]*MovedParticle.ori[2]);
                MovedParticle.ori[0] *= norm;
                MovedParticle.ori[1] *= norm;
                MovedParticle.ori[2] *= norm;
            }

            for( int vv = 0; vv < Config.part.size() && !inside; vv++ ){
                if(randP==vv) continue;
                if(vv < Config.Nc ){
                    if(overlapSP(MovedParticle,Config.part[vv], Config.l)) inside=true; 
                }else{
                    if(overlapSPSPH(MovedParticle,Config.part[vv], Config.l)) inside=true;
                }
            }

        }else{
            MovedParticle.pos[0] = MovedParticle.pos[0] + (uni(gen)-0.5)*pos_lambda;
            MovedParticle.pos[1] = MovedParticle.pos[1] + (uni(gen)-0.5)*pos_lambda;
            MovedParticle.pos[2] = MovedParticle.pos[2] + (uni(gen)-0.5)*pos_lambda;

            for( int vv = 0; vv < Config.Nc && !inside; vv++ ){
                if(overlapSPSPH(Config.part[vv],MovedParticle, Config.l)) inside=true; 
            }

            if(!inside){
                sx = MovedParticle.pos[0]*Config.w[0];
                if(sx==Config.W[0]) sx = Config.W[0]-1;
                sy = MovedParticle.pos[1]*Config.w[1];
                if(sy==Config.W[1]) sy = Config.W[1]-1;
                sz = MovedParticle.pos[2]*Config.w[2];
                if(sz==Config.W[2]) sz = Config.W[2]-1;

                int sx_max = sx + 1;
                int sy_max = sy + 1;
                int sz_max = sz + 1;
                
                int sx_min = sx - 1;
                int sy_min = sy - 1;
                int sz_min = sz - 1;

                for (int ix=sx_min; ix <= sx_max && !inside; ix++){
                    int kx = ix;
                    if(kx == Config.W[0]) kx =0;
                    else if(kx == -1) kx = Config.W[0]-1;

                    for (int iy=sy_min; iy <= sy_max && !inside; iy++){
                        int ky = iy;
                        if(ky == Config.W[1]) ky =0;
                        else if(ky == -1) ky = Config.W[1]-1;
                        for (int iz=sz_min; iz <= sz_max && !inside; iz++){
                            int kz = iz; 
                            if(kz == Config.W[2]) kz =0;
                            else if(kz == -1) kz = Config.W[2]-1;

                            int k = kx + Config.W[0]*(ky + Config.W[1]*kz);
            
                            int vv = Config.head[k];
                            
                            while(vv != -1 && !inside){
                                if(randP != vv ){
                                    if(overlapSPH(MovedParticle,Config.part[vv], Config.l)) inside=true; 
                                }
                                vv = Config.link[vv];
                            }
                        }
                    }
                }
            }
        }



        if(!inside){
            acc++;

            if( MovedParticle.pos[0] > Config.l[0] ) MovedParticle.pos[0] -= Config.l[0];
            if( MovedParticle.pos[1] > Config.l[1] ) MovedParticle.pos[1] -= Config.l[1];
            if( MovedParticle.pos[2] > Config.l[2] ) MovedParticle.pos[2] -= Config.l[2];

            if( MovedParticle.pos[0] < 0 ) MovedParticle.pos[0] += Config.l[0];
            if( MovedParticle.pos[1] < 0 ) MovedParticle.pos[1] += Config.l[1];
            if( MovedParticle.pos[2] < 0 ) MovedParticle.pos[2] += Config.l[2];

            
            if( randP >= Config.Nc){
                int k = sx + Config.W[0]*(sy + Config.W[1]*sz); 

                if(MovedParticle.cell != k){
                    if(Config.head[MovedParticle.cell] != randP ){
                        int v = Config.head[MovedParticle.cell];
                        int vv = Config.link[v];
                        while(vv!=randP){ 
                            v = vv;
                            vv = Config.link[v];
                        }
                        
                        Config.link[v] = Config.link[randP];
                    }else{
                        Config.head[MovedParticle.cell] = Config.link[randP];
                    }
                    
                    Config.link[randP] = Config.head[k];
                    Config.head[k] = randP;
                    MovedParticle.cell = k;
                }
            }

            Config.part[randP] = MovedParticle;
        }
    }
}


void Compression_step(){	
    double Vn = log(Config.Vbox)*vproc;
    Vn = exp(Vn);
    double ln = Config.l[0]/pow(Vn,1.0/3.0);

    Config.w[0] *= ln;
    Config.w[1] *= ln;
    Config.w[2] *= ln;

    ln = 1/ln;

    Config.l[0] *= ln;
    Config.l[1] *= ln;
    Config.l[2] *= ln;


    for( int i=0; i < Config.part.size(); i++){
        for( int vv=0; vv < 3; vv++){ 
            Config.part[i].pos[vv] *= ln;
        }
    }
    bool inside = false;

    for( int v=0; v < Config.part.size()-1 && !inside; v++){
        for( int vv=v+1; vv < Config.part.size() && !inside; vv++){
            if(v < Config.Nc){
                if(vv < Config.Nc ){
                    if(overlapSP(Config.part[v],Config.part[vv], Config.l)) inside=true; 
                }else{
                    if(overlapSPSPH(Config.part[v],Config.part[vv], Config.l)) inside=true; 
                }
            }else{
                if(vv < Config.Nc ){
                    if(overlapSPSPH(Config.part[vv],Config.part[v], Config.l)) inside=true; 
                }else{
                    if(overlapSPH(Config.part[v],Config.part[vv], Config.l)) inside=true; 
                }
            }
        }
    }

    if(!inside){ 
        Config.Vbox = Vn;
        Config.rhoN = N/Config.Vbox; 
        Config.rhoV = Config.rhoN*Config.Vsys; 
        vproc = 0.5*(3*vproc-1);
        if(vproc < 0.9 ) vproc=0.9;
        std::cout << "Compession successful! New rho = " << Config.rhoV << " and new vproc=" << vproc << std::endl;
        std::cout << "(lx,ly,lz) = (" << Config.l[0] << "," << Config.l[1] << "," << Config.l[2] << ")"<< std::endl;

        Config.write("Save/Config.dat",1);
        Config.write("Config.dat",0);
    }
    else{
        ln = 1/ln;
        Config.l[0] *= ln;
        Config.l[1] *= ln;
        Config.l[2] *= ln;

        for( int i=0; i < Config.part.size(); i++){
            for( int vv=0; vv < 3; vv++){ 
                Config.part[i].pos[vv] *= ln;
            }
        }

        RenewList();

        vproc = 0.5*(1+vproc);
        if(vproc > 0.9999){ 
            vproc = 0.9999;
            initCompress = true;
        }
        std::cout << "Compession unsuccessful! New vproc=" << vproc << std::endl;
        std::cout << "(lx,ly,lz) = (" << Config.l[0] << "," << Config.l[1] << "," << Config.l[2] << ")"<< std::endl;
        if (acceptance > 0.7){
            if( pos_lambda < maxpos) pos_lambda += 0.05;
            else pos_lambda = maxpos;
            if( ori_lambda < maxpos) ori_lambda += 0.05;
            else ori_lambda = maxpos;
            std::cout << "and changed pos_lambda=" << pos_lambda << " and ori_lambda=" <<  ori_lambda << std::endl;
        }
        if (acceptance < 0.5){
            if( pos_lambda < 0.01 ) pos_lambda *= 0.5;
            else pos_lambda -= 0.01;
            if( ori_lambda < 0.01 ) ori_lambda *= 0.5;
            else ori_lambda -= 0.01;
            std::cout << "and changed pos_lambda=" << pos_lambda << " and ori_lambda=" <<  ori_lambda << std::endl;
        }
    }
}
